#!/usr/bin/env python3
from sklearn.manifold import TSNE
import numpy as np
import pandas as pd
import torch
import seaborn as sns
from models_nn import load_mean_tensors
import matplotlib.pyplot as plt

np.random.seed(69)#la trik

train_data  = pd.read_csv("./train_tot.csv")
tensors_train = load_mean_tensors(train_data["seq_id"], "./train_tot__protein_esm2")
tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=10000)
tsne_results = tsne.fit_transform(tensors_train)

train_data['tsne-2d-one'] = tsne_results[:,0]
train_data['tsne-2d-two'] = tsne_results[:,1]
train_data.to_csv("train_tot_tsne.csv")

plt.figure(figsize=(16,10))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    hue="group",
    palette=sns.color_palette("hls", 10),
    data=train_data,
    legend="full",
    alpha=0.3
)
plt.show()
