#!/usr/bin/env python3


import numpy as np
import pandas as pd
import keras
from keras import Model
from keras.layers import Concatenate, Input, Dense
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize

import torch

def load_mean_tensors(prot_id_list, direct):
    res = []
    for i in range(len(prot_id_list)):
        cur_tensor = torch.load(direct + "/protein_id: " + str(prot_id_list[i]) + ".pt")["mean_representations"][33]
        res.append(cur_tensor.tolist())
    res_np = np.array(res)
    return res_np



# Data manipulation ------------------------------------------------------------
# dataframes
train_data = pd.read_csv("./train_all.csv")
test_data = pd.read_csv("./test.csv")
wildtype_data = pd.read_csv("./wildtypes.csv")

#embeddings
train_embedding = load_mean_tensors(train_data["seqid"], "./train_protein_esm2")
test_embedding  = load_mean_tensors(test_data["seq_id"], "./test_protein_esm2")
wildtype_embedding = load_mean_tensors(wildtype_data["gid"], "./wildtype_protein_esm2")

# Wildtype multiplication
# group_counts = train_data["gid"].value_counts().sort_index()

wt_mult = np.empty((train_data.shape[0], 1280), float)
for index, row in train_data.iterrows() :
    group_index_in_wt = int(row["gid"][1:])
    wt_mult[index, :] = wildtype_embedding[group_index_in_wt, :]


Y = np.array(train_data["dtm"])


# Model ------------------------------------------------------------------------

# simpler model
#
Y_normalized = normalize([Y])
X_train, X_test, Y_train, Y_test = train_test_split(train_embedding, Y_normalized, test_size=0.2)

inputs_spl = Input(shape=(1280,), name="embedding")
x_spl = Dense(32, "relu")(inputs_spl)
x_spl = Dense(16, "relu")(x_spl)
output_spl = Dense(1, "sigmoid")(x_spl)

model_spl = Model(inputs_spl, output_spl)
model_spl.compile(optimizer=keras.optimizers.Adam(learning_rate=0.1), loss="mean_squared_error", metrics=["mean_squared_error"])
history = model_spl.fit(X_train, Y_train, validation_split=0.2, epochs=50, batch_size=100)

def plot_loss(history):
  plt.plot(history.history['loss'], label='loss')
  plt.plot(history.history['val_loss'], label='val_loss')
  plt.xlabel('Epoch')
  plt.ylabel('Error [MPG]')
  plt.legend()
  plt.grid(True)

plot_loss(history)
plt.show()

predictions = model_spl.predict(X_test)
mse = np.mean((Y_test - predictions)**2)



def plot_modelspl(x, y):
  plt.scatter(x, y)
  plt.xlabel('Horsepower')
  plt.ylabel('MPG')
  plt.legend()
plot_modelspl(Y_test, predictions)
plt.show()

# Module muté

inputs_mt = Input(shape=1280)
x_mt = Dense(128, "relu")(inputs_mt)
out_mt = Dense(64, "relu")(x_mt)

# Module wt
inputs_wt = Input(shape=1280)
x_wt = Dense(128, "relu")(inputs_wt)
out_wt = Dense(64, "relu")(x_wt)

# modele tot
x = Concatenate()([x_mt, x_wt])
x = Dense(32, "relu")(x)
output = Dense(1, "sigmoid")(x)

model = Model([inputs_mt, inputs_wt], output)
model.compile(optimizer="sgd", loss="mse", metrics=["accuracy"])
history = model.fit([train_embedding, wt_mult], Y, validation_split=0.2, epochs=20, batch_size=50)

