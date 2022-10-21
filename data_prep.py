import Levenshtein
import matplotlib.pyplot as plt
import pandas as pd

train_data = pd.read_csv("./train.csv")

# filter data
# by size
train_data = train_data[(train_data["protein_sequence"].str.len() < 10000)
                        & (train_data["protein_sequence"].str.len() > 50)]
# by pH
train_data = train_data[(train_data["pH"] > 5) & (train_data["pH"] < 9)]

train_data = train_data.sort_index().reset_index(drop=True)

# make groups of proteins based on their pairwise Levenshtein distance
train_data_grouped = train_data
train_data_grouped["group"] = 0
last_group = 1
for i in range(train_data_grouped.shape[0]):
    for j in range(i+1, train_data_grouped.shape[0]):
        sequences = train_data_grouped.loc[(i, j), "protein_sequence"].tolist()
        if Levenshtein.distance(*sequences) <= 5:
            seq_group = list(
                filter(None, train_data_grouped.loc[(i, j), "group"]))
            if seq_group:
                group = min(seq_group)
            else:
                group = last_group
                last_group += 1
            train_data_grouped.loc[(i, j), "group"] = group

# select groups with 4 or more proteins
train_data_grouped = train_data.groupby("group").filter(lambda x: len(x) >= 4)
