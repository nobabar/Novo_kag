#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn as sk
from sklearn.linear_model import LinearRegression
import collections
train_data = pd.read_csv("./train.csv")
length_prots = [len(sequence) for sequence in train_data["protein_sequence"]]
length_prots = np.array(length_prots)
length_prots_reshaped = length_prots.reshape(-1, 1)

plt.hist(set(length_prots))
plt.show()
max_ind = np.argmax(length_prots)

min_ind = np.argmin(length_prots)
train_data.loc[min_ind, ]

# smol lr
reg = LinearRegression().fit(length_prots_reshaped, train_data["tm"])
reg.score(length_prots_reshaped, train_data["tm"])
inds_to = np.where(length_prots > 30000)
plt.scatter(np.delete(length_prots, inds_to), np.delete(np.array(train_data["pH"]), inds_to), alpha=0.01)
plt.show()

# phs
ph_idiot = np.where(train_data["pH"] > 14)
lengts_prots_ph_idiot = length_prots[ph_idiot[0]]
set(train_data.loc[ph_idiot[0],"data_source"])

#tm _distrib
plt.hist(train_data["tm"])
plt.show()

# doi
dois = collections.Counter(train_data["data_source"])

# mh_tm
plt.scatter(np.delete(np.array(train_data["pH"]), ph_idiot[0]), np.delete(np.array(train_data["tm"]), ph_idiot[0]), alpha=0.01)
plt.show()
