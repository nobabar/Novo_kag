#!/usr/bin/env python3


import itertools
import numpy as np
import pandas as pd
import keras
from keras import Input, Model
from keras.models import Sequential
from keras.layers import Activation, Add, Activation, BatchNormalization, AveragePooling1D
from keras.layers import Dense
from keras.layers import Dropout
from keras.layers import Conv1D, Flatten
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

import torch

def load_mean_tensors(prot_id_list, direct):
    res = []
    for i in range(len(prot_id_list)):
        cur_tensor = torch.load(direct + "/protein_id: " + str(prot_id_list[i]) + ".pt")["mean_representations"][33]
        res.append(cur_tensor.tolist())
        # tensors_series.append(cur_tensor)
    res_np = np.matrix(res)
    return res_np

def create_conc_XY(prot_data_frame, tensor_list):
    conc_X = []
    conc_Y = []
    groups = set(prot_data_frame["gid"])
    for group in groups:
        conc_prot = prot_data_frame[prot_data_frame["gid"] == group]
        for index in list(itertools.combinations(conc_prot.index,2)):
            tensor_conc = np.concatenate((
                np.squeeze(np.asarray(tensor_list[index[0]])),
                np.squeeze(np.asarray(tensor_list[index[1]]))
            ))
            diff_tm = conc_prot.loc[index[0], "dtm"] - conc_prot.loc[index[1], "dtm"]
            conc_Y.append(diff_tm)
            conc_X.append(tensor_conc)
    return (np.matrix(conc_X), np.array(conc_Y))


test_data = pd.read_csv("./test.csv")
train_data  = pd.read_csv("./train.csv")
thermo_data = pd.read_csv("./thermomut_grouped.csv", index_col=0)
max_seq_id = max(test_data["seq_id"])
thermo_index = list(range(max_seq_id + 1, len(thermo_data) + max_seq_id + 1))
thermo_data["seq_id"] = thermo_index

thermo_tensors = load_mean_tensors(thermo_data["seq_id"], "./thermo_protein_esm2")
conc_thermo_tensors, Y = create_conc_XY(thermo_data, thermo_tensors)


# Model


def dense_model():
    model = Sequential()
    model.add(Dense(32, input_dim=2560, activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(64, activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(1))
    return model

model_dense = dense_model()
model_dense.compile(optimizer="adam", loss="mse", metrics=["accuracy"], weighted_metrics=["accuracy"])
history_dense = model_dense.fit(conc_thermo_tensors, Y, validation_split = 0.2, epochs=20, batch_size=50)


plt.plot(history_dense.history['accuracy'])
plt.plot(history_dense.history['val_accuracy'])
plt.title('Model accuracy')
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='upper left')
plt.plot()
