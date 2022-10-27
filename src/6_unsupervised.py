#!/usr/bin/env python3
from sklearn.manifold import TSNE
import numpy as np
import pandas as pd
import torch
import seaborn as sns
from models_nn import load_mean_tensors
import matplotlib.pyplot as plt

np.random.seed(69)#la trik

train_data  = pd.read_csv("./train_all.csv", index_col=0)
wildtype_data = pd.read_csv("./wildtypes.csv")


tensors_train = load_mean_tensors(train_data["seqid"], "./train_protein_esm2")
tensors_wildtype = load_mean_tensors(wildtype_data["gid"], "./wildtype_protein_esm2")
tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=10000)
tsne_results = tsne.fit_transform(np.concatenate((tensors_train, tensors_wildtype)))

train_data['tsne_x'] = tsne_results[:,0]
train_data['tsne_y'] = tsne_results[:,1]
train_data.to_csv("train_tot_tsne.csv")
