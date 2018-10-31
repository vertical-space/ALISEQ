# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 15:17:59 2018

@author: cow082
"""

import sys

def gtf2bed(gtf, bed=None):
    """Convert a *.gtf file into a *.bed file
    
    Note: discard everything except exons
    """
    
    outlines = []
    count = 0
    
    with open(gtf, 'r') as v:
        for i, line in enumerate(v.readlines()):
            
            if line[0] == '#': # headers
                pass
                #print(line.strip())
                
            else:
                count += 1
                elements = line.strip().split('\t')
                
                if not len(elements) == 9:
                    #print(repr(elements))
                    pass
                
                elif elements[2] not in ['CDS', 'gene', 'exon', 'cDNA_match', 
                                         'lnc_RNA', 'mRNA', 'tRNA', 
                                         'origin_of_replication', 'pseudogene', 
                                         'region', 'transcript', 'V_gene_segment', 
                                         'rRNA', 'D_loop', 'snoRNA', 
                                         'snRNA', 'C_gene_segment', 'guide_RNA']:
                    #print(elements[2])
                    pass
                
                if elements[2] == 'exon':
                    meta = dict([item.split('=') for item in elements[8].split(';')])
                    data = [
                            elements[0],           # chromosome
                            elements[3],           # start
                            elements[4],           # end
                            elements[2],           # name ('exon')
                            meta.get('gene', '.'), # "score" (HGNC gene id)
                            elements[6],           # strand (+/-)
                            ]
                    outlines.append('\t'.join(data))
    
    if not bed:
        bed = gtf.split('.gtf')[0]+'.bed'
        
    with open(bed, 'w') as fOut:
        fOut.write('\n'.join(outlines))


def vcf2bed(vcf, bed=None, MIN_TOTAL_COUNT=10, MIN_ALLELE_COUNT=1):
    """Convert a *.vcf file output from bcftools::call into a *.bed file
    
    MIN_TOTAL_COUNT specifies the minimum number of high quality reads
    for a candidate to be retained
    
    MIN_ALLELE_COUNT sets a lower limit on the read depth of the 
    minor allele
    """
    
    outlines = []
    count = 0
    
    with open(vcf, 'r') as v:
        for i, line in enumerate(v.readlines()):
            
            if line[0] == '#': # headers
                #print(line.strip())
                pass
                
            else:
                count += 1
                elements = line.strip().split('\t')
                
                if not len(elements) == 10:
                    #print(repr(elements))
                    pass
                
                else:
                    metadata = dict([item.split('=') for item in elements[7].strip().split(';')])
                    #print(metadata)
                    """ 
                    # example:
                        
                    {'DP': '3', 'VDB': '0.66', 'SGB': '-0.453602', 'MQSB': '1', 
                    'MQ0F': '0', 'AF1': '1', 'AC1': '2', 'DP4': '0,0,1,1', 
                    'MQ': '60', 'FQ': '-32.988'}   
                    
                    DP4:    "Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases"

                    """
                    counts = metadata.get('DP4', '').split(',')
                    
                    if len(counts) == 4:
                        counts = [int(s) for s in counts]
                    else:
                        counts = [0,0,0,0]
                        
                    if sum(counts) < MIN_TOTAL_COUNT:
                        continue
                    
                    refCount = sum(counts[:2])
                    altCount = sum(counts[2:])
                                             
                    if min(refCount, altCount) < MIN_ALLELE_COUNT:
                        continue
                    
                    data = [
                            elements[0],             # chrom
                            elements[1],             # chromStart
                            str(int(elements[1])+1), # chromEnd
                            elements[3],             # refAllele			
                            elements[4],             # altAllele
                            str(refCount),           # refCount
                            str(altCount),           # altCount
                            ]
                    
                    outlines.append('\t'.join(data))
    
    if not bed:
        bed = vcf.split('.vcf')[0]+'.bed'
        
    with open(bed, 'w') as fOut:
        fOut.write('\n'.join(outlines))
        
#    print('\n'.join(outlines))


def csv2bed(csv, bed=None):
    """Convert a *.csv file output from allelecounter.py into a *.bed file
    """
    outlines = []
    headers = ['contig', 'position', 'variantID', 'refAllele', 'altAllele', 
               'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 
               'lowBaseQDepth', 'rawDepth', 'otherBases' ]
    
    with open(csv, 'r') as v:
        for i, line in enumerate(v.readlines()):
            
            if i == 0: # header
                cols = line.strip().split()
                if cols != headers:
                    print('\n'.join([
                            "Incorrect file type:",
                            "This function is solely designed to work on", 
                            "output files produced by allelecounter.py",
                            ]))
                    
                    sys.exit()
                
            else:
                elements = line.strip().split('\t')
                
                if not len(elements) == 12:
                    #print(repr(elements))
                    pass
                
                else:
                    data = [
                            elements[0],             # chrom
                            elements[1],             # chromStart
                            str(int(elements[1])+1), # chromEnd
                            elements[3],             # refAllele			
                            elements[4],             # altAllele
                            elements[5],             # refCount
                            elements[6],             # altCount
                            ]
                    outlines.append('\t'.join(data))
    
    if not bed:
        bed = csv.split('.csv')[0]+'.bed'
        
    with open(bed, 'w') as fOut:
        fOut.write('\n'.join(outlines))
    
#    with open(bed, 'r') as fIn:
#        print(fIn.read())

def main():

#    gtf2bed("C:/Users/cow082/water_buffalo_genome.gtf")
#    vcf2bed("C:/Users/cow082/ERR2353209.vcf")
#    csv2bed("C:/Users/cow082/ERR2353209.m2.csv")

    if sys.argv[1][-3:] == 'gtf':
        gtf2bed(sys.argv[1])
    elif sys.argv[1][-3:] == 'vcf':
        vcf2bed(sys.argv[1])
    elif sys.argv[1][-3:] == 'csv':
        csv2bed(sys.argv[1])

if __name__ == '__main__':
    main()

