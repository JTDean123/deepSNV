#!/usr/bin/env python3

# ----------------------------------
#
# Pipeline to train a CNN on SNVs
#
# Jason Dean, 12/17/2017
#
# ----------------------------------

import pandas as pd
import sys
import numpy as np
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from optparse import OptionParser
import feature_extract
import snvCNN


def getData(sample_num, vcf_path, genome_path, bam_path, len_path):
    '''
    Convert indexed, sorted bam file and VCF into CNN compatible data and labels

    Returns a 2-channel pileup/ref numpy array and a binary label array

    Parameters
    ----------
    sample_num : int
        number of snvs to sample for training
    vcf_path : str
        path to truth vcf file
    genome_path : str
        path to indexed reference genome
    bam_path : str
        path to indexed, sorted bam file
    len_path : str
        path to tab separated chromosome\tlength file


    Returns
    -------
    (X, y) : tuple
        X is a numpy array with dimensions [sample_num, 2, 5, 21]
        y is a numpy array of dimensions [sample_num * 2]
    '''

    # convert the vcf to a pandas df
    snps = feature_extract.vcf_to_pd(vcf_path)

    # randomly sample sample_num SNVs
    sys.stderr.write('INFO:  sampling {0} SNVs\n'.format(sample_num))
    snps_sample = snps.sample(sample_num)
    snps_sample = snps_sample.reset_index()

    # generate 2D matrix for snps
    snp_matrix = feature_extract.bam_matrix(bam_path, genome_path, snps_sample, 'snps')

    # sample non-snps
    # get the chromosome lengths
    chrom_len = feature_extract.get_genome(len_path)
    # randomly sample sample_num NON-SNVs
    nosnps = feature_extract.random_refs(sample_num, chrom_len, snps)

    # generate 2D matrix for non-snps
    nosnps_matrix = feature_extract.bam_matrix(bam_path, genome_path, nosnps, 'refs')

    # create feature vector.  SNV=1, nonSNV=0
    X = np.vstack((snp_matrix, nosnps_matrix))
    y = np.zeros(sample_num * 2)
    y[0:sample_num] = 1

    # save the feature and label data
    # TODO: make filename an arg
    sys.stderr.write('Saving feature and label data as snv.features.labels.npz\n')
    np.savez_compressed('snv.features.labels.npz', X, y)

    return (X, y)


def main():
    # read user passed arguments
    parser = OptionParser()
    parser.add_option("-n", "--sample_num", dest="sample_num", help="number of snps for training", default="100000", type='int')
    parser.add_option("-v", "--vcf_path", dest="vcf_path", help="path to vcf")
    parser.add_option("-g", "--genome_path", dest="genome_path", help="path to reference genome")
    parser.add_option("-b", "--bam_path", dest="bam_path", help="path to sorted, indexed bam file")
    parser.add_option("-l", "--len_path", dest="len_path", help="path to file containing chrom lengths")
    parser.add_option("-e", "--epochs", dest="epochs", help="number of training epochs, default 50", default='50', type='int')

    (options, args) = parser.parse_args()

    X, y = getData(options.sample_num, options.vcf_path, options.genome_path, options.bam_path, options.len_path)

    # split into test train sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

    # build a model
    trained_cnn = snvCNN.train_cnn(X_train, y_train, options.epochs)

    # evaluate on test data
    predictions = snvCNN.test_cnn(trained_cnn, X_test)

    # generate confusion matrix
    print('saving confusion matrix as deepSNV.confusion.csv')
    preds = pd.DataFrame(confusion_matrix(y_test, predictions))
    preds.columns = ['Predicted Non-SNV', 'Predicted SNV']
    preds.index = ['Actual Non-SNV', 'Actual SNV']
    preds.to_csv('deepSNV.confusion.csv')

    # evaluate accuracy on test data
    acc = snvCNN.accuracy(y_test, predictions)
    print('Accuracy on test data:  ', acc)

    # save the model
    print('Saving model as deepSNV.h5')
    trained_cnn.save('deepSNV.h5')

    print('Exiting')
    sys.exit(0)


# go time
if __name__ == "__main__":
    main()
