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
    # convert indexed, sorted bam file and VCF into CNN compatible format

    # sample snps
    snps = feature_extract.vcfSample(vcf_path)
    print('sampling', sample_num, 'SNPS')
    snps_sample = snps.sample(sample_num)
    snps_sample = snps_sample.reset_index()

    # generate 2D matrix for snps
    sample_type = 'snps'
    snp_matrix = feature_extract.bamMatrix(bam_path, genome_path, snps_sample, sample_type)

    # sample non-snps
    chrom_len = feature_extract.getGenome(len_path)
    nosnps = feature_extract.randomRefs(sample_num, chrom_len, snps)

    # generate 2D matrix for non-snps
    sample_type = 'refs'
    nosnps_matrix = feature_extract.bamMatrix(bam_path, genome_path, nosnps, sample_type)

    # create feature vector.  SNV=1, nonSNV=0
    X = np.vstack((snp_matrix, nosnps_matrix))
    y = np.zeros(sample_num * 2)
    y[0:sample_num] = 1

    # save the feature and label data
    print('Saving feature and label data as snv.features.labels.npz')
    np.savez_compressed('snv.features.labels.npz', X, y)

    return (X, y)


def main():

    # read user passed arguments
    parser = OptionParser()
    parser.add_option("-n", "--sample_num", dest="sample_num", help="number of snps for training", default="10000")
    parser.add_option("-v", "--vcf_path", dest="vcf_path", help="path to vcf")
    parser.add_option("-g", "--genome_path", dest="genome_path", help="path to reference genome")
    parser.add_option("-b", "--bam_path", dest="bam_path", help="path to sorted, indexed bam file")
    parser.add_option("-l", "--len_path", dest="len_path", help="path to file containing chrom lengths")
    parser.add_option("-e", "--epochs", dest="epochs", help="number of training epochs, default 50", default='50')

    (options, args) = parser.parse_args()

    sample_num = int(options.sample_num)
    vcf_path = options.vcf_path
    genome_path = options.genome_path
    bam_path = options.bam_path
    len_path = options.len_path
    epochs = int(options.epochs)

    snp_data = getData(sample_num, vcf_path, genome_path, bam_path, len_path)
    X, y = snp_data

    # split into test train sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

    # build a model
    trained_cnn = snvCNN.trainCNN(X_train, y_train, epochs)

    # evaluate on test data
    predictions = snvCNN.testCNN(trained_cnn, X_test)

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