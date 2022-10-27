#!/usr/bin/env python3

from numpy import index_exp
import pandas as pd


# Save data as fasta ---------------------------
def save_data_as_fasta(data, output_file, prot_seq_name, seqid = None):
    with open(output_file, "w") as fileout:
        for ind, prot in data.iterrows():
            if seqid:
                name = prot[seqid]
            else:
                name = ind
            fileout.write(f'>protein_id: {name}\n')
            fileout.write(f'{prot[prot_seq_name]}\n')


# Using train all

test_data = pd.read_csv("./test.csv")
save_data_as_fasta(test_data, "test.fasta", "seq_id")
train_data  = pd.read_csv("./train_all.csv", index_col=0)
wildtype_data = pd.read_csv("./wildtypes.csv")
save_data_as_fasta(train_data, "train.fasta", "protein_sequence", "seqid")
save_data_as_fasta(wildtype_data, "wildtype.fasta", "wildtype", "gid")
