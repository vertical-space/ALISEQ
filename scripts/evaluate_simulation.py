# -*- coding: utf-8 -*-

"""

.. Created on Mon Oct  8 16:21:06 2018

   @author: cow082

   evaluate_simulation.py

"""

import argparse
import pandas as pd


def evaluate(args):
    '''examine a collection of results.csv files
    '''
    stats = pd.read_csv(args.statsFile)
    
    #print(stats.columns)
    """
    Index(['Unnamed: 0', 'gene_id', 'transcript_snp_id', 'ase', 'rpkm',
       'total_read_count', 'ref_read_count', 'alt_read_count',
       'transcript_length', 'relative_snp_position', 'ref_nucleotide',
       'alt_nucleotide', 'average_coverage', 'snps_per_transcript'],
      dtype='object')

    """
    
    total_true_positives = set(stats[                                         \
                (stats['ase'] == True) &                                      \
                (stats['total_read_count'] >= 1) &                            \
                (stats['snps_per_transcript'] >= 1)                           \
                ].gene_id)
    
    total_true_negatives = set(stats[                                         \
                (stats['ase'] == False) &                                     \
                (stats['total_read_count'] >= 1) &                            \
                (stats['snps_per_transcript'] >= 1)                           \
                ].gene_id)
    
    for f in args.infileList.split(','):
        
        results = {}
        
        hits = set(pd.read_csv(f, sep='\t').gene)
        
        multiplier = int(f.split('_x')[1].split('.summary.csv')[0])
        #e.g. results/simulation_7_x16.summary.csv   --> 16
        
        valid_true_positives = set(stats[                                     \
                    (stats['ase'] == True) &                                  \
                    (stats['total_read_count']*multiplier >= args.minreads) & \
                    (stats['snps_per_transcript'] >= args.minsnps)            \
                    ].gene_id)
        
        valid_true_negatives = set(stats[                                     \
                    (stats['ase'] == False) &                                 \
                    (stats['total_read_count']*multiplier >= args.minreads) & \
                    (stats['snps_per_transcript'] >= args.minsnps)            \
                    ].gene_id)
    
        results['TP'] = total_true_positives.intersection(hits)
        results['FP'] = total_true_negatives.intersection(hits)
        results['TN'] = total_true_negatives.difference(hits)
        results['FN'] = total_true_positives.difference(hits)
        
        results['TPE'] = valid_true_positives.intersection(hits)
        results['FPE'] = valid_true_negatives.intersection(hits)
        results['TNE'] = valid_true_negatives.difference(hits)
        results['FNE'] = valid_true_positives.difference(hits)
        
        print('simulation_x%d' % multiplier)
        print('\n'.join(['   %s %d' % (k, len(v)) for k,v in results.items()]))
        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infileList', required=True,
                        help="specify input filenames as a comma-seperated list")
    parser.add_argument('-s', '--statsFile', required=True,
                        help="file containing sample statistics")
    parser.add_argument('-r', '--minreads', required=False, default=10, type=int,
                        help="minimum average coverage for a gene to be considered")
    parser.add_argument('-n', '--minsnps', required=False, default=5, type=int,
                        help="minimum number of SNPs for a gene to be considered")
    evaluate(parser.parse_args())


if __name__ == '__main__':
    main()

