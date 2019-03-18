# ----------------------------------
#
# Predict SVN status via trained CNN
#
# Jason Dean, 12/17/2017
#
# ----------------------------------

import pandas as pd
from optparse import OptionParser
import sys
import feature_extract
import keras
import snvCNN


def locations(loc_path):

    # open locations for prediction
    try:
        locs = pd.read_csv(loc_path, index_col=False, header=None)
        locs.columns = ['chrom', 'pos']
        return locs
    except:
        sys.stderr('ERROR opening location file, exiting')
        sys.exit(1)


def loadModel(model_path):

    # open trained model for prediction
    try:
        cnn = keras.models.load_model(model_path)
        return cnn
    except:
        sys.stderr('ERROR opening trained model, exiting')
        sys.exit(1)


def writeOutput(predictions, coords):

    # write output
    out_preds = []
    for i in predictions:
        if i == 1:
            out_preds.append('YES')
        else:
            out_preds.append('NO')
    coords['Predicted SNV'] = out_preds

    try:
        coords.to_csv('predictions.csv')
        print('Successfully wrote output to predictions.csv')
    except:
        sys.stderr('ERROR writing output, exiting')
        sys.exit(1)


def main():

    # read user passed arguments
    parser = OptionParser()
    parser.add_option("-p", "--preds_path", dest="preds_path", help="path to csv with chromosome and coordinate for SNV prediction")
    parser.add_option("-g", "--genome_path", dest="genome_path", help="path to reference genome")
    parser.add_option("-b", "--bam_path", dest="bam_path", help="path to sorted, indexed bam file")
    parser.add_option("-m", "--model_path", dest="model_path", help="path to trained CNN model")


    (options, args) = parser.parse_args()
    preds_path = options.preds_path
    genome_path = options.genome_path
    bam_path = options.bam_path
    model_path = options.model_path

    # load coordinates for prediction
    coordinates = locations(preds_path)

    # load trained model
    model = loadModel(model_path)

    # generate feature matrix
    sample_type = 'prediction'
    pred_matrix = feature_extract.bam_matrix(bam_path, genome_path, coordinates, sample_type)

    # make predictions
    print('Making predictions....')
    preds = snvCNN.test_cnn(model, pred_matrix)

    # write output
    writeOutput(preds, coordinates)

    # goodbye
    print('Exiting')


# go time
if __name__ == "__main__":
    main()