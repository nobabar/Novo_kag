#!/usr/bin/env python3
from sklearn.manifold import TSNE
import numpy as np
import pandas as pd
import torch
import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(69)  # la trik


def load_mean_tensors(prot_id_list, direct):
    res = []
    for i in range(len(prot_id_list)):
        cur_tensor = torch.load(
            direct + "/protein_id_ " + str(prot_id_list[i]) + ".pt")["mean_representations"][33]
        res.append(cur_tensor.tolist())
    res_np = np.array(res)
    return res_np


train_data = pd.read_csv("./data/train_all.csv", index_col=0)
wildtype_data = pd.read_csv("./data/wildtypes.csv")
test_data = pd.read_csv("./data/test.csv")


tensors_train = load_mean_tensors(
    train_data["seqid"], "./data/train_protein_esm2")
tensors_wildtype = load_mean_tensors(
    wildtype_data["gid"], "./data/wildtype_protein_esm2")
tensors_test = load_mean_tensors(
    test_data["seqid"], "./data/test_protein_esm2")
tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=10000)
tsne_results = tsne.fit_transform(np.concatenate(
    (tensors_train, tensors_wildtype, tensors_test)))


train_data['tsne_x'] = tsne_results[:len(train_data), 0]
train_data['tsne_y'] = tsne_results[:len(train_data), 1]

wildtype_data['tsne_x'] = tsne_results[len(
    train_data):len(train_data)+len(wildtype_data), 0]
wildtype_data['tsne_y'] = tsne_results[len(
    train_data):len(train_data)+len(wildtype_data), 1]

test_data['tsne_x'] = tsne_results[len(train_data)+len(wildtype_data):, 0]
test_data['tsne_y'] = tsne_results[len(train_data)+len(wildtype_data):, 1]
train_data.to_csv("./data/train_tot_tsne.csv")
wildtype_data.to_csv("./data/wildtype_tsne.csv")
test_data.to_csv("./data/test_tsne.csv")
