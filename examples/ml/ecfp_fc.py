#! /usr/bin/env python

import numpy as np
from keras.models import Sequential
from keras.optimizers import Adam
from keras.layers.core import Dense
from rmgpy.cnn_framework.data import prepare_full_train_data_from_multiple_datasets
import logging

def lr(epoch):
	lr = 0.0007
	return float(lr* np.exp(- epoch / 30.0))

# hyper-params
hidden_units = 50
momentum = 0.5
input_dim = 512
epochs = 150

# Load the datasets
datasets =  [
			('rmg', 'sdata134k', 'polycyclic_2954_table', 0.1), 
			('rmg', 'sdata134k', 'cyclic_O_only_table', 0.1)
			]
X_test, y_test, X_train, y_train = prepare_full_train_data_from_multiple_datasets(datasets, using_ecfp=True)

#################################
## Model specification

## Start from an empty sequential model where we can stack layers
model = Sequential()

model.add(Dense(hidden_units, input_dim=input_dim, activation='tanh'))

model.add(Dense(1, activation='linear'))

##################################

## Compile the model with categorical_crossentrotry as the loss, and stochastic gradient descent (learning rate=0.001, momentum=0.5,as the optimizer)
model.compile(loss='mse', optimizer=Adam(lr=0.01))

## Fit the model (10% of training data used as validation set)
best_val_loss = 9e10
no_improve_count = 0
for epoch in range(epochs):
	model.optimizer.lr.set_value(lr(epoch))

	print('\nEpoch {}/{}, lr = {}'.format(epoch+1, epochs, lr(epoch)))

	history_callback = model.fit(np.array(X_train), 
								np.array(y_train), 
								nb_epoch=1, 
								batch_size=1,
								validation_split=0.1)
	loss_history = history_callback.history
	loss = loss_history['loss'][0]
	val_loss = loss_history['val_loss'][0]
	if val_loss < best_val_loss:
		best_val_loss = val_loss
		no_improve_count = 0
	else:
		no_improve_count += 1
		print('{} epochs without inner_val_loss progress'.format(no_improve_count))

	if no_improve_count >= 10:
		break

## Evaluate the model on test data
objective_score = model.evaluate(np.array(X_test), np.array(y_test), batch_size=32)

# objective_score is a tuple containing the loss as well as the accuracy
print ("\nLoss on test set:"  + str(objective_score))

model.save_weights('ecfp_fc.h5', overwrite = True)