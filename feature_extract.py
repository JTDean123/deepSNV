# ----------------------------------
#
# Functions to convert a sorted, indexed bam file into
# a format suitable for CNN model training.  SNPs are
# determined from a user provided VCF.
#
# Jason Dean, 12/17/2017
#
# ----------------------------------


import pysam
import pandas as pd
from pysam import VariantFile
import sys
import numpy as np
import random


def vcf_to_pd(path):
    """
    Convert a VCF file to a pandas df

    Return a pandas df with cols [chrom, pos, ref, alts]

    Parameters
    ----------
    path : str
        path to VCF file containing ground truth mutation calls

    Returns
    -------
    vcfs : pandas df
        vcfs is a dataframe with dimensions vcf_variant calls x 4

    """
    # function to randomly sample a vcf file

    vcfs = pd.DataFrame()

    try:
        print('loading SNP VCF')
        vcf = VariantFile(path)
    except:
        sys.stderr('ERROR:  Unable to load VCF')
        sys.exit(1)

    chrom = []; pos = []; ref = []; alts = []

    for rec in vcf.fetch():
        pos.append(rec.pos)
        chrom.append(rec.chrom)
        ref.append(rec.ref)
        alts.append(str(rec.alts).split(',')[0].split('(')[-1])

    vcfs['chrom'] = chrom
    vcfs['pos'] = pos
    vcfs['ref'] = ref
    vcfs['alts'] = alts

    return vcfs


def bam_matrix(bam_path, ref_path, coords, sample_type):
    """
    Create a [coords.shape[0], 2, 5, 22] feature matrix

    Parameters
    ----------
    bam_path : str
        path to sorted, indexed bam file
    ref_path : str
        path to indexed reference file
    coords : pandas df
        pandas df - [chrom, pos, ref, alt]
    sample_type : str
        either 'snps' or 'refs', only used for printing status updates

    Returns
    -------
    snps : numpy array
        a numpy array with dimensiuons [coords.shape[0], 2, 5, 22]

    """

    sys.stderr.write('INFO:  creating feature matrix for {0}\n'.format(sample_type))
    observations = coords.shape[0]
    snps = np.zeros((observations, 2, 5, 22))

    # nucleotide lookup dictionary
    nt_dict = {'a': 0,
               't': 1,
               'c': 2,
               'g': 3,
               'n': 4}

    # load bam file
    try:
        sys.stderr.write('INFO: loading BAM file\n')
        bam = pysam.Samfile(bam_path, 'rb')
    except:
        sys.stderr.write('ERROR:  Unable to load BAM\n')
        sys.exit(1)

    # load reference genome
    try:
        sys.stderr.write('INFO: loading reference genome\n')
        ref = pysam.Fastafile(ref_path)
    except:
        sys.stderr.write('ERROR:  Unable to load reference genome\n')
        sys.exit(1)

    sys.stderr.write('INFO:  Creating nucleotide feature matrix\n')
    # get a nucleotide specific matrix
    for i, j in enumerate(coords['pos']):
        for k, coordinate in enumerate(list(range(j - 11, j + 11))):

            # fill in reference matrix
            nt_ref = ref.fetch(coords.loc[i, 'chrom'].split('chr')[-1], coordinate, coordinate + 1).lower()
            if nt_ref == '': nt_ref = 'n'
            snps[i, 0, nt_dict[nt_ref], k] += 1

            for read in bam.pileup(coords.loc[i, 'chrom'], j - 1, j):
                if read.pos == int(coordinate):
                    for pileupread in read.pileups:
                        # fill in sequence read matrix
                        if pileupread.query_position is not None:
                            nt_snp = pileupread.alignment.query_sequence[pileupread.query_position].lower()
                        else:
                            nt_snp = 'n'
                        snps[i, 1, nt_dict[nt_snp], k] += 1

    sys.stderr.write('INFO:  generated feature matrix for {0}\n'.format(sample_type))

    return snps


def get_genome(path):
    """
    Get the chromosome lengths for the reference genome

    TODO:  could extract this from the fai ref file

    Parameters
    ----------
    path : str
        path to chromosome lengths file.  chr1\tchr1_length\nchr2\tchr_2_length ect

    Returns
    -------
    chrom ref : pandas df
        pandas data frame with dimensions

    """
    # get ref genome length
    chrom19 = []
    len19 = []
    try:
        with open(path) as file:
            line = file.read().splitlines()
            for i in line:
                chrom19.append(i.split('\t')[0])
                len19.append(int(i.split('\t')[1]))
    except:
        sys.stderr('ERROR:  Unable to load reference genome lengths')
        sys.exit(1)

    chrom_ref = pd.DataFrame()
    chrom_ref['chrom'] = chrom19
    chrom_ref['length'] = len19

    return chrom_ref


def random_refs(num_pick, chrom_len, pdvcf):
    """
    Pick random non-snv spots in the genome

    Returns a pandas df containing chrom and reference nt seq

    Parameters
    ----------
    num_pick : int
        number of coordinates to sample
    chrom_len : pandas df
        pandas df with columns [chrom, length]
    pdvcf : pandas df
        pandas df of vcf, created by vcfSample()

    Returns
    -------
    no_snp : pandas df
        pandas df with dimensions [num_pick, 2] and cols [chrom, no_snps]

    """

    # the following two blocks will generate a df with three columns: chrom, snv counts, and fraction
    # calculate the number of SNVs on each chromosome
    snps_counts = pdvcf.groupby('chrom')['pos'].nunique()
    chrom38 = list(chrom_len['chrom'])
    snps_counts = snps_counts[chrom38]
    snps_counts.columns = ['Chromosome', 'Number of SNVs']

    # calculate fraction of SNVs in each genome
    chroms = snps_counts.index.tolist()
    counts = snps_counts.iloc[:, ]
    counts = list(counts.fillna(0))
    svns_distribution = pd.DataFrame()
    svns_distribution['chrom'] = chroms
    svns_distribution['counts'] = counts
    svns_distribution['fraction'] = svns_distribution['counts'] / svns_distribution['counts'].sum()

    # determine the amount of ref coords to select from each chromosome
    svns_distribution['samples'] = np.random.multinomial(num_pick,
                                                         list(svns_distribution['fraction']),
                                                         size=1).flatten()

    # generate a pseudo VCF dataframe with random ref coordinates that are not SNVs
    sys.stderr.write('INFO:  generating random reference non-snp locations\n')
    no_snps = pd.DataFrame()
    no_snps_chrom = []
    no_snps_spot = []
    padding = 10000

    # TODO: this is too hard to read.
    # TODO: better way is to subset pdvfc by chrom and check 'random_spot is in pdvcf['pos']'
    # iterate through each chromosome
    for index, spot in enumerate(svns_distribution['samples']):
        if spot == 0: continue

        # get length of the current chromosome
        chrom = svns_distribution['chrom'][index]
        length = chrom_len[chrom_len['chrom'] == chrom]['length']
        # don't sample from the first and last 10kb b/c these are N in ref genome
        region = list(range(0 + padding, int(length) - padding))

        for sample in range(spot):
            random_spot = random.choice(region)
            is_a_snp = True

            while True:
                # is the randomly selected coordinate a SNV? if not, keep it.  if yes, break
                if pdvcf[(pdvcf['chrom'] == chrom) & (pdvcf['pos'] == random_spot)].shape[0] == 1:
                    random_spot = random.choice(region)
                else:
                    break

            no_snps_chrom.append(chrom)
            no_snps_spot.append(random_spot)

    no_snps['chrom'] = no_snps_chrom
    no_snps['pos'] = no_snps_spot

    return no_snps
