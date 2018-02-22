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


def trainCNN(X_train, y_train, epochs):
    # reshape the data
    print('Training model with:', X_train.shape[0], 'training observations and ', epochs, 'epochs')
    X_train = X_train.reshape(X_train.shape[0], 5, 22, 2)

    # flex the model
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

    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=['accuracy'])
    model.fit(X_train, y_train, batch_size=16, epochs=epochs, verbose=1)

    return model


def testCNN(model, X_test):
    # evaluate the model on the test data
    print('Evaluating model with:', X_test.shape[0], 'testing observations')
    X_test = X_test.reshape(X_test.shape[0], 5, 22, 2)
    predictions = model.predict_classes(X_test)

    return predictions


def accuracy(y_test, predict):
    # calcualate total accuracy
    right = 0
    for i, j in enumerate(predict):
        if j == y_test[i]:
            right += 1
    accuracy = right / len(y_test)

    return accuracy
