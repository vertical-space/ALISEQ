#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

.. Created on Thu Sep 20 13:53:33 2018

   @author: chris

"""

import os
import argparse
import logging
import random
import numpy as np
import scipy.stats as stats
import pandas as pd
from randomdict import RandomDict
from Bio import Seq


"""

.. Note: the first 3 functions to generate a random sequence library, 
   are not presently used, but could be useful later

"""

def createRandomSequence(mu=3392, sigma=2600, m=300, M=10000):
    '''create a random ATCG string with random length. 
    Note: use the normal distribution to choose the length (mean: 3392bp, s.d 2600bp), 
    but discard lengths outside the range [m,M]
    '''
    length = 0
    while not m <= length <= M:
        length = int(np.random.normal(mu, sigma))
    return ''.join([random.choice('ATCG') for i in range(length)])


def createRandomSequenceIndex(n):
    '''constructs a dictionary {seq_id: sequence} of random ATCG sequences 
    (for default values see the randomGene() function
    '''
    return dict([('gene_%s' % (i+1), createRandomSequence()) for i in range(n)])


def createRandomSequenceLibrary(n_seqs=400):
    '''generate a random sequence library and save to file
    '''
    seqs = createRandomSequenceIndex(n_seqs)
    transcriptome = fasta(seqs)
    with open('simulated_transcriptome.fasta', 'w') as F:
        F.write(transcriptome)


def fastWrap(text, width=None):
    '''much quicker than textwrap.wrap
    '''
    if not width:
        return text
    
    return '\n'.join(text[i:i+width] for i in range(0, len(text), width))


def fasta(seq_index, wrap=None):
    '''format a dictionary of {label:sequence} pairs into a fasta string
    '''
    return '\n'.join(['>%s\n%s' % (k, fastWrap(v, wrap)) for k,v in seq_index.items()])


def fastq(seq_index, args, wrap=None):
    '''format a dictionary of {label:sequence} pairs into a fastq string
    '''
    return '\n'.join(['@%s\n%s\n+\n%s' % (k, fastWrap(v, wrap), \
                fastWrap(qualScore(len(v), args), wrap)) for k,v in seq_index.items()])


def qualScore(i, args):
    '''select a quality string at random from a real RNA-Seq dataset
    
    NOTE: It is assumed that desired strings are <= the length of the 
    available strings.
    
    There is an assert statement in createReadLibrary() which enforces this 
    at the time of args.qualdata construction.
    '''
    return args.qual_scores.random_value()[:args.read_length+1]
    
    
def randomDecision(probability=0.5):
    '''returns a random True/False based on:
            <probability>       chance of returning True
            1 - <probability>   chance of returning False
    '''
    assert 0 <= probability <= 1
    
    return random.random() < probability


def parseHeaders(s):
    '''custom fasta header parser specific to this dataset
    
    might not work with other datasets
    '''
    parts = s.split('gene=')
    if len(parts) > 1: 
        return parts[1].split(']')[0]
    else:
        return s.split('[')[0].strip()


def indexFasta(fasta, limit=None):
    '''re-written to account for some entries having '>' in the title string!
    '''
    index = {}
    with open(fasta) as F:
        head = ''
        tail = []
        records = 0
        for line in F.readlines():
            if line[0] == '>':
                if tail:
                    if records >= limit:
                        break
                    seq = ''.join(tail)
                    if head in index:
                        if len(seq) > len(index[head]): # keep only the longest transcript for each gene
                            index[head] = seq
                            records += 1
                    else:
                        index[head] = seq
                        records += 1
                head = parseHeaders(line.lstrip('>').strip())
                tail = []
            else:
                tail.append(line.strip())
    return index


def createAltAllele(seq, SNP_frequency):
    '''simulate the mutation process to generate a heterozygous allele
    '''
    
    ref_seq = list(seq)
    
    # duplicate the refseq
    alt_seq = ref_seq.copy()
    
    # choose randomised positions to mutate
    mutation_indices = mutateString(ref_seq, SNP_frequency)
    
    # get the ref nucleotides at these positions
    ref_nucleotides = [ref_seq[i] for i in mutation_indices]
    
    # choose random nucleotides for each position. make sure they are different from the ref nucleotides
    alt_nucleotides = [random.choice(list(set(list('ATCG')).difference(set(list(r))))) for r in ref_nucleotides]
    
    #perform the substitutions in the alt_seq
    for index, ref, alt in zip(mutation_indices, ref_nucleotides, alt_nucleotides):
        alt_seq[index] = alt
    
    # return the mutated sequence plus the details of any alt/ref bases it contains
    return ''.join(alt_seq), zip(mutation_indices, ref_nucleotides, alt_nucleotides)


def mutateString(ref_seq, SNP_frequency):
    '''based on a specified frequency, return a list of positions to mutate
    
    Note: there are many possible ways this could be done
    '''
    i = len(ref_seq)
    t = 0
    mutation_indices = []
    while t < i:
        p = np.random.geometric(SNP_frequency)
        t += p
        if t < i:
            mutation_indices.append(t)
    return mutation_indices


def prepareData(args):
    '''index the reference transcriptome, assign read counts, then decide 
    which genes are going to be ASE
    '''

    assert args.output_tag # blank tags not allowed
    assert os.path.exists(args.rpkm_file)
    assert os.path.exists(args.transcriptome)
    assert os.path.exists(args.fastq)

    # index the reference expression data
    RPKM_index = pd.read_csv(args.rpkm_file)
    n_seqs = RPKM_index.shape[0] # number of rows in the df
    
    # extract sequence data from an existing reference library, e.g. a transcriptome file
    seqs = indexFasta(args.transcriptome, limit=args.gene_limit)
    heads, tails = list(zip(*seqs.items()))
    n_seqs = min([n_seqs, len(seqs)])
    
    # randomly obtain a list of RPKM values with a realistic distribution profile
    rpkm = random.sample(list(RPKM_index[args.tissue]), n_seqs)
    
    '''alternative models for testing'''
#    rpkm = [10]*n_seqs
#    rpkm = sorted(list(RPKM_index[args.tissue]), reverse=True)[:n_seqs]
    
    # randomly select ASE genes at the desired frequency
    decisions = [randomDecision(args.ase_frequency) for i in range(n_seqs)]

    df = pd.DataFrame({
        'gene':heads,
        'seq' :tails,
        'ase' :decisions,
        'rpkm':rpkm,
        })

    # determine the total read count per locus from the RPKM values, desired library size and gene lengths
    df['reads'] = df.apply(lambda row: \
        round(row.rpkm * (len(row.seq) / 1000) * (args.lib_size / 1000000)), 
        axis=1)
    
    #print(np.sum(df.reads))
    
    return df


def defineCountsAndSequences(df, args):
    '''loop through each gene (transcript) constructing a second
    allele (based on a defined mutation rate) and read counts for each allele 
    (divide up the total reads between ref & alt alleles: proportions depend on whether the 
    gene has been assigned as ASE or not)
    '''
    statistics = {}
    snp_counter = 1
    ratio = 1/args.inbalance, 1 # e.g. (0.5, 1)
    proportions = [i/sum(ratio) for i in ratio] # e.g. (0.33, 0.66)
    
    # make a copy of the dataframe
    df1 = df.copy()
    
    # add some empty columns to the new dataframe
    df1['alt_seq'] = ''
    df1['ref_read_count'] = 0
    df1['alt_read_count'] = 0

    # iterate over each gene/transcript
    for row in df.itertuples():
        
        # create the alt allele
        alt_seq, stats = createAltAllele(row.seq, args.snp_frequency)
        
        # determine allele-specific read counts
        if row.ase:
            ref_read_count, alt_read_count = [round(row.reads*p) for p in proportions]
        else:
            ref_read_count, alt_read_count = [round(row.reads*p) for p in [0.5, 0.5]]
        
        # add the allele data to the new dataframe    
        df1.loc[row[0], 'alt_seq'] = str(alt_seq)
        df1.loc[row[0], 'ref_read_count'] = int(ref_read_count)
        df1.loc[row[0], 'alt_read_count'] = int(alt_read_count)
        
        genewise_snp_count = 1 # only used to assign <transcript_snp_id>
        
        # format stats for logging
        for pos, ref, alt in stats:
            statistics['SNP_%s' % snp_counter] = pd.Series({
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
    
    # process and format the stats output
    statistics = pd.DataFrame(statistics).T
    statistics['average_coverage'] = statistics.total_read_count * args.read_length / statistics.transcript_length
    snp_count = dict([(gene, list(statistics.gene_id).count(gene)) for gene in list(set(statistics.gene_id))])
    statistics['snps_per_transcript'] = statistics.apply(lambda row: snp_count[row.gene_id], axis=1)

    return df1, statistics


def chooseInsert(seq, args, READ_COUNTER, max_attempts=5):
    '''choose the insert start and end positions based on a 
    predetermined distribution of insert sizes
    '''
    start = random.choice(range(len(seq)))
    
    if args.uniform_insert_sizes: # try to make all reads the same length
        
        insert_size = min([
                args.insert_mean, # this is our target size
                len(seq) - start, # but reads sometimes need to be smaller if they are close to the end
                ])
                
    else: # we are going to use a (potentially skewed) normal distribution to determine read sizes
        
        '''make sure reads do not exceed the specified max length, 
        but also that they do not exceed the end of the transcript
        '''
        insert_max = min([args.insert_max, len(seq)-start]) 
        
        '''if the maximum insert size is smaller than the smallest value in the
        insert_distribution, it is not possible to use random sampling to get a 
        suitable value, so use one we already have'''
        if insert_max <= args.insert_distribution_min:
            insert_size = insert_max
            
        else:
            count = 0
            while True:
                insert_size = args.insert_distribution.random_value()
                if args.insert_min <= insert_size <= insert_max:
                    break
                '''Avoid spending too much time sampling. If the insert_max 
                is close to the margins of the insert_distribution, it may take 
                millions of attempts to locate a suitable value!'''
                if count == max_attempts: # give up after n failed attempts!
                    insert_size = insert_max
                    break
                count += 1
                
    end = int(start+insert_size)
    
    return start, end
            

def invert(i, s):
    '''return a coordinate corresponding to a position on the reverse strand
    that matches a coordinate on the forward strand, i.e.
    
    0     5
    | --> | 
    AAAAAAAAAAAAAAAA
    TTTTTTTTTTTTTTTT
          |        |
          9    <-- 0

    '''
    return len(s) - i - 1


#def pad(s, args):
#    '''if a read is shorter than intended, pad with N's
#    '''
#    return s + 'N'*max(0, args.read_length-len(s))


def createReadPair(seq, reads_file_1, reads_file_2, READ_COUNTER, args):
    '''simulate a single pair of reads in an NGS experiment
    '''   
    
    # select an insert
    start, end = chooseInsert(seq, args, READ_COUNTER) # forward strand position
    insert = seq[start:end+1]
    template_1 = insert + args.adaptor + ''.join('N'*args.read_length)
    read_1 = template_1[:args.read_length]
    template_2 = str(Seq.Seq(insert).reverse_complement()) + args.adaptor + ''.join('N'*args.read_length)
    read_2 = template_2[:args.read_length]
    
    # ensure the first read starts on line 0
    if READ_COUNTER == 0: 
        linebreak = ''
    else:
        linebreak = '\n' 

    reads_file_1.write(linebreak+fastq({
        'simulation '+str(READ_COUNTER)+'/1': read_1,
        }, args))
    reads_file_2.write(linebreak+fastq({
        'simulation '+str(READ_COUNTER)+'/2': read_2,
        }, args))
    
    if READ_COUNTER % 1000 == 0: # write buffer to disk every 1000 reads
        reads_file_1.flush()
        reads_file_2.flush()
        
    READ_COUNTER += 1
    
    return READ_COUNTER
    

def readQualscores(args):
    '''basic function for extracting quality strings from a fastq file
    '''
    quals = []
    
    with open(args.fastq, 'r') as fastq:
        for i, line in enumerate(fastq.readlines()):
            if i % 4 == 3: # process the 3rd line out of every 4 lines 
                line = line.strip()
                if len(line) >= args.read_length:
                    quals.append(line)
                    
    random.shuffle(quals)
    
    return quals


def createReadLibrary(df, reads_file_1, reads_file_2, args):
    '''loop through a dataframe containing sequences (ref and alt alleles) 
    and expression levels (ref and alt counts) and create two dictionaries of 
    simulated NGS reads (fwd and reverse read pairs)
    '''
    READ_COUNTER = 0
    args.read_target = np.sum(df.ref_read_count) + np.sum(df.alt_read_count)
    
    insert_distribution = [int(i) for i in list(stats.skewnorm.rvs(
            a=args.skew,
            size=args.read_target+10000, # add a buffer in case of rounding errors etc.
            loc=args.insert_mean,
            scale=args.insert_sd,
            ))]
    
    random.shuffle(insert_distribution)
    
    args.insert_distribution = RandomDict([(k,v) for k,v in enumerate(insert_distribution)])
    
    args.insert_distribution_min = min(insert_distribution)
    
    args.qual_scores = RandomDict([(k,v) for k,v in enumerate(readQualscores(args))])
                
    # iterate over each gene/transcript
    for row in df.itertuples():
                        
        # create ref_seq reads
        for i in range(row.ref_read_count):
            READ_COUNTER = createReadPair(
                row.seq, reads_file_1, reads_file_2, READ_COUNTER, args)

        # create alt_seq reads
        for i in range(row.alt_read_count):
            READ_COUNTER = createReadPair(
                row.seq, reads_file_1, reads_file_2, READ_COUNTER, args)
            

def simulateSequencingExperiment(args):
    '''create a simulated NGS dataset containing allelic inbalance in the 
    expression of a defined proportion of genes
    
    this function can be modified to output multiple sequencing simulations 
    from a single simulated RNA sample, but with varying read counts to 
    assess the effect of increasing sequencing depth
    '''

    logging.basicConfig(filename=args.logfile, level=logging.DEBUG)
    for arg, value in vars(args).items():
        logging.info('%s: %r' % (arg, value))

    # simulate generation of the RNA library
    df = prepareData(args)
    df, statistics = defineCountsAndSequences(df, args)
    outfile_1 = 'simulation_%s.stats.txt'  % args.output_tag
    if not args.clobber: 
        assert not os.path.exists(outfile_1)
    statistics.to_csv(outfile_1)
    
    # simulate the sequencing run/s
    multipliers = [1,]
    if args.use_multipliers:
        multipliers.extend([2,4,8,16,32])
    for m in multipliers:
        df_new = df.copy()
        df_new.ref_read_count *= m
        df_new.alt_read_count *= m
        outfile_2 = 'simulation_%s_x%d.1.fastq' % (args.output_tag, m)
        outfile_3 = 'simulation_%s_x%d.2.fastq' % (args.output_tag, m)
        for f in [outfile_2, outfile_3]:
            if not args.clobber:
                assert not os.path.exists(f)
        reads_file_1 = open(outfile_2, 'w')
        reads_file_2 = open(outfile_3, 'w')
        createReadLibrary(df_new, reads_file_1, reads_file_2, args)
        reads_file_1.close()
        reads_file_2.close()


def createTestInsert(s_fwd, args):
    '''create an insert (for testing purposes only)
    '''
    start_1 = random.choice(range(len(s_fwd) - args.insert_min + 1))
    start_2 = chooseInsert(s_fwd, start_1, args) # forward strand position
    return start_2 - start_1


def printStats(a, args):
    print('\t'.join([str(obj) for obj in ['', 'target', 'actual']]))
    print('\t'.join([str(obj) for obj in ['mean', args.insert_mean, np.mean(a)]]))
    print('\t'.join([str(obj) for obj in ['median', '-', np.median(a)]]))
    print('\t'.join([str(obj) for obj in ['sd', args.insert_sd, np.std(a)]]))
    print('\t'.join([str(obj) for obj in ['min', args.insert_min, np.min(a)]]))
    print('\t'.join([str(obj) for obj in ['max', args.insert_max, np.max(a)]]))
    
    
def test_1(args):
    '''plot insert sizes from an invariant pseudo-transcriptome
    '''
    
    import matplotlib.pyplot as plt

    a = np.array([createTestInsert('A'*10000, args) for i in range(5000)])
    plt.hist(a, bins=50)
    plt.show()
    printStats(a, args)
    

def test_2(args):
    '''plot insert sizes from an actual transcriptome
    '''
    
    import matplotlib.pyplot as plt
    
    seqs = indexFasta(args.transcriptome, limit=args.gene_limit)
    heads, tails = list(zip(*seqs.items()))
    a = np.array([createTestInsert(random.choice(tails), args) for i in range(5000)])
    plt.hist(a, bins=50)
    plt.show()
    printStats(a, args)
    

def test_3(args):
    '''
    model the insert size distribution of an actual transcriptome
    
    >>> from scipy import stats

    # choose some parameters
    >>> a, loc, scale = 1.3, -0.1, 2.2
    # draw a sample
    >>> sample = stats.skewnorm(a, loc, scale).rvs(1000)
    
    # estimate parameters from sample
    >>> ae, loce, scalee = stats.skewnorm.fit(sample)
    >>> ae
    1.2495366661560348
    >>> loce
    -0.039775813819310835
    >>> scalee
    2.1126121580965536
    
    NOTE: the file 'ERR2353209.insert_size_metrics.txt' has insert data based on 
    mapping RNA-Seq reads to a genome. Therefore, some reported inserts will 
    span splice junctions and thus be recorded as longer than they actually are. 
    It was generated by running picard::CollectInsertSizeMetrics
    
    I do not plan to modify it, however one way would be to map reads to the 
    transcriptome instead, then re-measure the insert sizes.
    '''
    
    import matplotlib.pyplot as plt
    
    df = pd.read_csv('ERR2353209.insert_size_metrics.txt', 
                     delim_whitespace=True, header=0, skiprows=10)
    
    sample = []
    for row in df.itertuples():
        for i in range(int(row[2])):
            sample.append(row[1])
        
    ae, loce, scalee = stats.skewnorm.fit(sample)
    # print(ae, loce, scalee)
    # 8.4 97.3 149.2
    # set defaults to: 8, 100, 150
        
    insert_distribution = stats.skewnorm.rvs(a=args.skew, 
                                             loc=args.insert_mean, 
                                             scale=args.insert_sd,
                                             size=args.lib_size,
                                             )
    plt.hist(insert_distribution, bins=50)
    plt.show()
    printStats(insert_distribution, args)
    
    
    insert_distribution = stats.skewnorm.rvs(a=ae,
                                             loc=loce, 
                                             scale=scalee,
                                             size=args.lib_size,
                                             )
    plt.hist(insert_distribution, bins=50)
    plt.show()
    printStats(insert_distribution, args)
    
    
def test_4(args):
    '''
    '''
    simulateSequencingExperiment(args)
    

def argParser():
    
    parser = argparse.ArgumentParser()
    
    # required argument
    parser.add_argument('-a', '--output_tag', required=True,
                        help='a tag that gets added to all output file names')
    # optional args...
    parser.add_argument('-b', '--snp_frequency', required=False, type=float, default=0.000625, 
                        help='proportion of coding bases that are heterozygous')
    parser.add_argument('-c', '--ase_frequency', required=False, type=float, default=0.1, 
                        help='proportion of genes with allelic expression inbalance')
    parser.add_argument('-d', '--lib_size', required=False, type=int, default=2315672, 
                        help='target library size (in read pairs)')
    parser.add_argument('-e', '--read_length', required=False, type=int, default=124, 
                        help='length of single reads')
    parser.add_argument('-f', '--inbalance', required=False, type=float, default=2,
                        help='extent of skew toward most abundant allele')
    parser.add_argument('-g', '--gene_limit', required=False, type=int, default=1000000,
                        help='only process the first <gene_limit> transcripts')
    parser.add_argument('-i', '--adaptor', required=False, default='AGATCGGAAGAGC',
                        help='e.g. illumina 3p adaptor')
    parser.add_argument('-j', '--rpkm_file', required=False, default='rpkm_fixed.csv',
                        help='a table of RPKM expression values')
    parser.add_argument('-k', '--tissue', required=False, default='testis',
                        help='which tissue to model expression values on?')
    parser.add_argument('-l', '--transcriptome', required=False, default="water_buffalo_transcriptome.fasta",
                        help='a transcriptome file containing the longest transcript at each protein coding locus')
    parser.add_argument('-m', '--use_multipliers', action='store_true',
                        help='whether to simulate multiple runs with different read depths')
    parser.add_argument('-n', '--clobber', action='store_true',
                        help='whether to overwrite existing files')
    parser.add_argument('-o', '--insert_mean', required=False, type=float, default=100, 
                        help='average insert length')
    parser.add_argument('-p', '--insert_sd', required=False, type=float, default=150, 
                        help='standard deviation of insert length')
    parser.add_argument('-q', '--skew', required=False, type=float, default=8, 
                        help='skews the insert size distribution')
    parser.add_argument('-r', '--insert_min', required=False, type=int, default=0, 
                        help='minimum allowed insert size')
    parser.add_argument('-s', '--insert_max', required=False, type=int, default=1000, 
                        help='maximum allowed insert size')
    parser.add_argument('-t', '--uniform_insert_sizes', action='store_true',
                        help='make all inserts equal to <insert_mean> size')    
    parser.add_argument('-u', '--logfile', required=False, default='logfile.txt',
                        help='make all inserts equal to <insert_mean> size')
    parser.add_argument('-v', '--fastq', required=False, default='ERR2353209.1.fastq',
                        help='an example fastq file that we will steal quality strings from')
    
    return parser

    
def main():

    myParser = argParser()
    
    simulateSequencingExperiment(myParser.parse_args())

#    random.seed(18)
#    np.random.seed(20)
#
#    test_1(myParser.parse_args([
#            '-a', 'test_1', 
#            '--insert_min', '125',
#            '--insert_max', '3000', 
#            '--insert_mean', '180',
#            '--insert_sd', '50',
##            '--uniform_insert_sizes', 
#            ]))

#    test_2(myParser.parse_args([
#            '--output_tag', 'test_2', 
#            '--transcriptome', 'C:/Users/cow082/2genes.fa',
#            '--insert_min', '0',
#            '--insert_max', '1000', 
#            '--insert_mean', '180',
#            '--insert_sd', '50',
##            '--uniform_insert_sizes', 
#            ]))
#    
#    test_3(myParser.parse_args([
#            '--output_tag', 'test_3', 
#            '--transcriptome', 'C:/Users/cow082/2genes.fa',
#            '--lib_size', '1000',
#            '--clobber',
#            '--insert_min', '0',
#            '--insert_max', '1000', 
#            '--insert_mean', '180',
#            '--insert_sd', '75',
#            '--skew', '4',
#            ]))
#    
#    test_4(myParser.parse_args([
#            '--output_tag', 'test_4', 
#            '--transcriptome', 'C:/Users/cow082/2genes.fa',
#            '--lib_size', '10000',
#            '--clobber',
#            ]))


if __name__ == '__main__':
    main()

