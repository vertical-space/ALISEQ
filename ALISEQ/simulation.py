#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 13:53:33 2018

@author: chris
"""

import random
#from collections import Counter
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import math
from Bio import Seq

random.seed(42)

def random_gene(m=300,M=10000):
    """create a random ATCG string with length randomly chosen from within the range [m,M]
    """
#    length = random.choice(range(m,M+1))
    length = 0
    while not m < length < (M+1):
        length = int(np.random.normal(3392,2600)) # mu, sigma
    return ''.join([random.choice('ATCG') for i in range(length)])


def fastwrap(text, width):
    """much quicker than textwrap.wrap"""
    return '\n'.join(text[i:i+width] for i in range(0, len(text), width))


def fasta(seq_index, wrap=80):
    """format a dictionary of {label:sequence} pairs into a fasta string
    """
    return '\n'.join(['>%s\n%s' % (k, fastwrap(v, wrap)) for k,v in seq_index.items()])


def fastq(seq_index, wrap=80):
    """format a dictionary of {label:sequence} pairs into a fasta string
    """
    return '\n'.join(['>%s\n%s\n+\n%s' % (k, fastwrap(v, wrap), fastwrap('F'*len(v), wrap)) for k,v in seq_index.items()])


def generate_random_seqs(n):
    return dict([('gene_%s' % (i+1), random_gene()) for i in range(n)])


def decision(probability):
    return random.random() < probability


def parse_headers(s):
    parts = s.split('gene=')
    if len(parts) > 1: 
        return parts[1].split(']')[0]
    else:
        return s.split('[')[0].strip()
    

def index_fasta(fasta, limit=None):
    """loads the whole file into memory. 
    could easily be implemented as a generator.
    """
    index = {}
    with open(fasta) as F:
        records = F.read().strip().split('>')[1:]
        for rec in records[:limit]:
            lines = rec.strip().split('\n')
            head = parse_headers(lines[0].strip())
            tail = ''.join([L.strip() for L in lines[1:]])
            if tail:
                index[head] = tail
    return index

def todo():
    SNP_frequency = 0.000625 # in an average human, 0.0625% of the exome is heterozygous (18798 SNPs / 30 MB)
    ase_frequency = 0.05 # 5% of genes have ASE
    library_size = 2315672 # 2.3M reads (target only - the actual number of reads will be different)
    read_length = 100 # 100bp PE reads
    adaptor_3p = 'GGGGTTTT'
    adaptor_5p = 'AAAACCCC'
    reads_file_1 = "simulated_reads_1.fastq"
    reads_file_2 = "simulated_reads_2.fastq"
    stats_file = "stats.csv"
    inbalance = 2 # alt alleles have 2x abundance compared to ref alleles
    ratio = 1/inbalance, 1 # e.g. (0.5, 1)
    proportions = [i/sum(ratio) for i in ratio] # e.g. (0.33, 0.66)
    
    """generate a random sequence library"""
    #n_seqs=400
    #seqs = generate_random_seqs(n_seqs)
    #transcriptome = fasta(seqs)
    #with open("simulated_transcriptome.fasta", 'w') as F:
    #    F.write(transcriptome)
        
    """use an existing reference sequence library, e.g. a transcriptome file"""
    input_fasta_file = "water_buffalo_transcriptome.fasta"
    seqs = index_fasta(input_fasta_file, limit=4000) # first 4000 transcripts only
    n_seqs = len(seqs)
    
    
    heads, tails = list(zip(*seqs.items()))
    decisions = [decision(ase_frequency) for i in range(n_seqs)] # randomly select ASE genes at the desired frequency
    #print(Counter(decisions))
    #Counter({False: 90, True: 10})
    
    RPKM_index = pd.read_csv('rpkm.csv')
    rpkm_ref = random.sample(list(RPKM_index['RPKM']), n_seqs) # randomly obtain a list of RPKM values with a realistic distribution profile
    
    count = 0
    reads_1 = {}
    reads_2 = {}

def make_read(s):
    global count
    
    start_1 = random.choice(range(len(s)))
    start_2 = random.choice(range(len(s)))
    
    str_s = ''.join(s)
    s_fwd = str_s + adaptor_3p + ''.join('N'*read_length)
    reads_1[str(count)+"/1"] = s_fwd[start_1:start_1+read_length]
    
    s_rev = str(Seq.Seq(str_s).reverse_complement())
    s_rev = s_rev + adaptor_5p + ''.join('N'*read_length)
    reads_2[str(count)+"/2"] = ''.join(s_rev[start_2:start_2+read_length])
    
    count += 1

def predict_SNP_count():
    """An example of weighted random choosing (currently unused)
    only works in py3.6
    """
    return random.choices(
            population=range(5),
            weights=[0.1, 0.1, 0.2, 0.3, 0.3],
            k=1
            )


def mutate(ref_seq):
    """once positions to mutate are chosen, identify the ref nucleotides at 
    these positions, then select replacements at random. Construct and return 
    the mutated sequence plus a list of records documenting the mutations: 
        [(position, ref_nucleotide, alt_nucleotide), ...]
    """
    alt_seq = ref_seq.copy()
#    mutations = random.choice(range(1,21)) # the number of SNPs is not uniformly distributed - this should be biased towards fewer SNPs
    #mutation_indices = random.sample(range(len(row.seq)), mutations)
    mutation_indices = mutate2(ref_seq)
    ref_nucleotides = [ref_seq[i] for i in mutation_indices]
    alt_nucleotides = [random.choice(list(set(list('ATCG')).difference(set(list(r))))) for r in ref_nucleotides]
    for index, ref, alt in zip(mutation_indices, ref_nucleotides, alt_nucleotides):
        alt_seq[index] = alt
    return alt_seq, zip(mutation_indices, ref_nucleotides, alt_nucleotides)


def mutate2(ref_seq):
    """randomly assign positions to mutate by using the negative binomial 
    distribution to scan L-R along th egene length skipping n nucleotides 
    between successive mutations until the end of the sequence is reached.
    """
    i = len(ref_seq)
    t = 0
    mutation_indices = []
    while t < i:
        p = np.random.negative_binomial(1, SNP_frequency)
        t += p
        if t < i:
            mutation_indices.append(t)
    return mutation_indices


def main():
    
    df = pd.DataFrame({
        'gene':heads,
        'seq':tails,
        'ase':decisions,
        'rpkm':rpkm_ref,
        })

    plt.plot(sorted(df.rpkm))
    
    plt.save('rpkm.png')
    #print(sorted([math.floor(row.rpkm) for row in df.itertuples() if row.ase]))
    
    """determine the total read count per locus from the RPKM values, desired library size and gene lengths"""
    df['reads'] = df.apply(lambda row: round(row.rpkm * (len(row.seq) / 1000) * (library_size / 1000000)), axis=1)
    
    print(np.sum(df.reads))
    
    mutation_statistics = {}
    snp_counter = 1
    
    for row in df.itertuples(): # iterate over each gene/transcript
        
        ref_seq = list(row.seq)
        alt_seq, stats = mutate(ref_seq)
        
        # determine allele-specific read counts
        if row.ase:
            ref_read_count, alt_read_count = [round(row.reads*p) for p in proportions]
        else:
            ref_read_count, alt_read_count = [round(row.reads*p) for p in [0.5, 0.5]]
        
        # write reads
        for r in range(ref_read_count):
            make_read(ref_seq)
        for r in range(alt_read_count):
            make_read(alt_seq)
            
        genewise_snp_count = 1
        for pos, ref, alt in stats:
            mutation_statistics['SNP_%s' % snp_counter] = pd.Series({
                    'gene_id' : row.gene, 
                    'transcript_snp_id': '%s_%d' % (row.gene, genewise_snp_count),
                    'ase': row.ase,
                    'rpkm': row.rpkm, 
                    'total_read_count': row.reads, 
                    'ref_read_count': ref_read_count, 
                    'alt_read_count': alt_read_count, 
                    'transcript_length': len(row.seq), 
                    'relative_snp_position': pos, 
                    'ref_nucleotide': ref, 
                    'alt_nucleotide': alt,
                    })
            snp_counter += 1
            genewise_snp_count += 1
    
    with open(reads_file_1, 'w') as F:
        F.write(fastq(reads_1))
    
    with open(reads_file_2, 'w') as F:
        F.write(fastq(reads_2))

    df = pd.DataFrame(mutation_statistics).T
    df['average_coverage'] = df.total_read_count * read_length / df.transcript_length
    snp_count = dict([(gene, list(df.gene_id).count(gene)) for gene in list(set(df.gene_id))])
    df['snps_per_transcript'] = df.apply(lambda row: snp_count[row.gene_id], axis=1)
    df.to_csv(stats_file)


if __name__ == "__main__":
    main()
