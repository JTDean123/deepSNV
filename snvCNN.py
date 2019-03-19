# ----------------------------------
#
# Build and train and test a CNN to predict if a
# genomic coordinate in a bam file is a SNP.
# Feature data should be generated via feature_extract.py
#
# Jason Dean, 12/17/2017
#
# ----------------------------------


from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers.convolutional import Convolution2D, MaxPooling2D
import sys


def train_cnn(X_train, y_train, epochs):
    """
    Train a CNN using 2-channel ngs derived pileup data

    return a trained model

    Parameters
    ----------
    X_train : numpy array
        sample feature data dimensions [len(y_train), 2, 5, 21]
    y_train : numpy array
        feature labels
    epochs : int
        number of epochs for training

    Returns
    -------
    model : h5
        keras (h5) hierarchical data format trained model

    """
    # reshape the data
    print('Training model with:', X_train.shape[0], 'training observations and ', epochs, 'epochs')
    X_train = X_train.reshape(X_train.shape[0], 5, 22, 2)

    # build and train the modl
    model = Sequential()

    model.add(Convolution2D(32, kernel_size=(2, 4),
                            strides=(1, 1),
                            activation='relu',
                            input_shape=X_train.shape[1:]))
    model.add(MaxPooling2D(pool_size=(1, 1), strides=(1, 1)))

    model.add(Convolution2D(64, kernel_size=(2, 4),
                            strides=(1, 1), activation='relu',
                            input_shape=X_train.shape[1:]))
    model.add(MaxPooling2D(pool_size=(2, 2), strides=(1, 1)))

    model.add(Flatten())

    model.add(Dropout(0.2))
    model.add(Dense(500, activation='relu'))

    model.add(Dropout(0.2))
    model.add(Dense(32, activation='relu'))

    model.add(Dropout(0.2))
    model.add(Dense(1, activation='sigmoid'))

    # TODO - add a callback to kill training if loss stops decreasing

    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=['accuracy'])
    model.fit(X_train, y_train, batch_size=16, epochs=epochs, verbose=1)

    return model


def test_cnn(model, X_test):
    """
    Make predictions on test data

    Parameters
    ----------
    model
    X_test : numpy array
        test data with dimensions [number_of_test_samples, 2, 5, 21]
    Returns
    -------
    predictions : numpy array
        a np array with dimensions [X_test.shape[0], 1]
        entries will either be 1 or 0

    """
    # evaluate the model on the test data
    sys.stderr.write('Evaluating model with: {0} test observations'.format(X_test.shape[0]))
    X_test = X_test.reshape(X_test.shape[0], 5, 22, 2)
    predictions = model.predict_classes(X_test)

    return predictions


def accuracy(y_test, predict):
    """
    Calculate accuracy

    Parameters
    ----------
    y_test : numpy array
        labels of test data
    predict : numpy array
        output of test_cnn

    Returns
    -------
    accuracy : float
        (number of correct predictions on test samples) / (number of test samples)

    """
    # calculate total accuracy
    right = 0
    for i, j in enumerate(predict):
        if j == y_test[i]:
            right += 1
    accuracy = right / len(y_test)

    return accuracy
