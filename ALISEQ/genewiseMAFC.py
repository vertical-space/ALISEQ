# -*- coding: utf-8 -*-

"""

.. Created on Wed Aug 29 12:37:31 2018

   @author: cow082

"""

import numpy as np
import pandas as pd
import math
import argparse
import sys


def meanAbsoluteFoldChange(frame):
    """find the average log2 fold change of SNPs (alleles) for a single gene
    """

    #print('meanAbsoluteFoldChange: run...', frame)
   
    frame['fc'] = frame.apply(
            lambda row: row.a_altCount / row.a_refCount, axis=1)
    
    frame['log2fc'] = frame.apply(
            lambda row: math.log(row.fc, 2), axis=1)
    
    frame['absLog2fc'] = frame.apply(
            lambda row: abs(row.log2fc), axis=1)
    
    frame['aveLog2fc'] = np.mean(frame['absLog2fc'])
    
    return pd.Series({
            'snp_count': int(len(frame)),
            'ave_read_depth': round((np.sum(frame.a_altCount) + np.sum(frame.a_refCount)) / len(frame), 1),
            'aveLog2fc': frame['aveLog2fc'].iloc[0],
            })

    
def enrichASE(args):
    """identify genes with differential expression between two alleles
    
    requires: annotated SNP data (*.merged.csv) output from bedtools intersect
    
    How it works:
        1. group SNPs together by gene
        2. For each SNP in that gene, calculate the absolute(log2(fold change))
        3. take the average of absolute(log2(fold change)) values across all SNPS for that gene
        4. filter on:
            - average read depth at variable positions for the gene
            - number of SNPs in the gene
            - average(absolute(log2(fold change))) for the gene
    
    Two optional args with no defaults change the program behaviour:
        1. --output <file> sends the results to file rather than stdout
        2. --gene <gene> yields data for every SNP in a specific gene. Without
           this arg, the program produces instead a summary table of all 
           significant genes in order of priority 
           (sorted by: average(absolute(log2(fold change)))) 
           
    """
    
    # the input has no header, so supply one
    headers = ['a_chrom', 'a_chromStart', 'a_chromEnd', 'a_refBase', 
               'a_altBase', 'a_refCount', 'a_altCount', 
               'b_chrom', 'b_chromStart', 'b_chromEnd', 'b_featureType', 
               'b_gene', 'b_orientation']
    
    # read the input data
    df1 = pd.read_csv(args.input, sep='\t', header=None, names=headers)
    
    # calculate meanAbsoluteFoldChange per gene
    df2 = df1.groupby(by='b_gene').apply(meanAbsoluteFoldChange)
    
    df2 = df2.sort_values(['aveLog2fc'], ascending=False)

    # filter hits
    df2 = df2[df2['ave_read_depth'] >= args.minCoverage]
    df2 = df2[df2['snp_count'] >= args.minSnps]
    df2 = df2[df2['aveLog2fc'] >= math.log(args.minFoldChange, 2)]
    
    df2.index.name = 'gene'
    
    # filter the main data table to identify expression at the per-SNP level 
    # for the significant genes identified above
    hits = pd.Series({'gene': sorted(df2.index)})
    df3 = df1[df1.b_gene.isin(hits.gene)]
    df3.set_index('b_gene', inplace=True)
    df3.index.name = 'gene'
    df3 = df3.drop(['b_chrom', 'b_chromStart', 'b_chromEnd', 
                    'b_featureType', 'b_orientation'], axis=1)
    df3.columns = ['chrom', 'start', 'end', 'refBase', 
                   'altBase', 'refCount', 'altCount']
    
    # choose where to send the output
    if args.output:
        save_location = args.output
    else:
        save_location = sys.stdout
        
    if args.gene:
        # print specific details for a single target gene
        df3.loc[args.gene,:].to_csv(save_location, sep='\t')
    else:
        # display a summary of all hits
        df2.head(len(df2)).to_csv(save_location, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, 
                        help="annotated SNP data from bedtools intersect")
    parser.add_argument('-r', '--minCoverage', required=False, type=int, default=10, 
                        help="minimun average read depth at variable sites")
    parser.add_argument('-s', '--minSnps', required=False, type=int, default=5, 
                        help="minimum number of SNPs in a given gene")
    parser.add_argument('-f', '--minFoldChange', required=False, type=float, default=2,
                        help="minimum average absolute fold change")
    parser.add_argument('-g', '--gene', required=False, 
                        help="gene of interest")
    parser.add_argument('-o', '--output', required=False, 
                        help="output file")
    enrichASE(parser.parse_args())

if __name__ == "__main__":
    main()

