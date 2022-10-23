#!/usr/bin/env python3

import pandas as pd


# Save data as fasta ---------------------------
def save_data_as_fasta(data, output_file):
    with open(output_file, "w") as fileout:
        for _, prot in data.iterrows():
            fileout.write(f'>protein_id: {prot["seq_id"]}\n')
            fileout.write(f'{prot["protein_sequence"]}\n')


test_data = pd.read_csv("./test.csv")
train_data  = pd.read_csv("./train.csv")

save_data_as_fasta(test_data, "test.fasta")
save_data_as_fasta(train_data, "train.fasta")

