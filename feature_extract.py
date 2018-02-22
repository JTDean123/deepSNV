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


def vcfSample(path):
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

    return (vcfs)


def bamMatrix(bam_path, ref_path, coords, sample_type):
    # convert a bam file into two 5x21x(#obs) matrices via provided coordinates

    print('creating feature matrix for', sample_type)
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
        print('loading BAM file')
        bam = pysam.Samfile(bam_path, 'rb')
    except:
        sys.stderr('ERROR:  Unable to load BAM')
        sys.exit(1)

    # load reference genome
    try:
        print('loading reference genome')
        ref = pysam.Fastafile(ref_path)
    except:
        sys.stderr('ERROR:  Unable to load reference genome')
        sys.exit(1)

    print('Creating nucleotide feature matrix')
    # get a nucleotide specific matrix
    for i, j in enumerate(coords['pos']):
        for k, coordinate in enumerate(list(range(j - 11, j + 11))):

            # fill in reference matrix
            nt_ref = ref.fetch(coords.loc[i, 'chrom'].split('chr')[-1], coordinate, coordinate + 1).lower()
            #nt_ref = ref.fetch(coords.loc[i, 'chrom'], coordinate, coordinate + 1).lower()
            if nt_ref == '': nt_ref = 'n'
            snps[i, 0, nt_dict[nt_ref], k] += 1

            for read in bam.pileup(coords.loc[i, 'chrom'], j - 1, j):
                if read.pos == int(coordinate):
                    for pileupread in read.pileups:
                        # print(pileupread.alignment.query_sequence[pileupread.query_position])
                        # fill in sequence read matrix
                        if pileupread.query_position is not None:
                            nt_snp = pileupread.alignment.query_sequence[pileupread.query_position].lower()
                        else:
                            nt_snp = 'n'
                        snps[i, 1, nt_dict[nt_snp], k] += 1

    print('generated feature matrix for', sample_type)

    return snps


def getGenome(path):
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


def randomRefs(num_pick, chrom_len, pdvcf):
    # get num_pick number of random (sort of) spots in the genome
    # calculate the number of snvs on each chromosome

    snps_counts = pd.DataFrame()
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

    # sample from a genome wide multinomial distribution
    svns_distribution['samples'] = np.random.multinomial(num_pick,
                                                         list(svns_distribution['fraction']),
                                                         size=1).flatten()

    # generate a VCF with random samples
    print('generating random reference non-snp locations')
    no_snps = pd.DataFrame()
    no_snps_chrom = []
    no_snps_spot = []
    padding = 10000

    for index, spot in enumerate(svns_distribution['samples']):
        if spot == 0: continue

        # get length of the current chromosome
        chrom = svns_distribution['chrom'][index]
        length = chrom_len[chrom_len['chrom'] == chrom]['length']
        region = list(range(0 + padding, int(length) - padding))

        for sample in range(spot):
            random_spot = random.choice(region)
            is_a_snp = True

            while True:
                if pdvcf[(pdvcf['chrom'] == chrom) & (pdvcf['pos'] == random_spot)].shape[0] == 1:
                    random_spot = random.choice(region)
                else:
                    break

            no_snps_chrom.append(chrom)
            no_snps_spot.append(random_spot)

    no_snps['chrom'] = no_snps_chrom
    no_snps['pos'] = no_snps_spot

    return no_snps