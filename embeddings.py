#!/usr/bin/env python3

from numpy import index_exp
import pandas as pd


# Save data as fasta ---------------------------
def save_data_as_fasta(data, output_file, seqid = None):
    with open(output_file, "w") as fileout:
        for ind, prot in data.iterrows():
            if seqid:
                name = seqid
            else:
                name = ind
            fileout.write(f'>protein_id: {name}\n')
            fileout.write(f'{prot["protein_sequence"]}\n')


# test_data = pd.read_csv("./test.csv")
# train_data  = pd.read_csv("./train_grouped.csv", index_col=0)
# thermo_data = pd.read_csv("./thermomut_grouped.csv", index_col=0)

# # Add thermo id
# max_seq_id = max(test_data["seq_id"])
# thermo_index = list(range(max_seq_id + 1, len(thermo_data) + max_seq_id + 1))
# thermo_data["seq_id"] = thermo_index


# # Add train id
# max_seq_id = max(thermo_data["seq_id"])
# train_index = list(range(max_seq_id + 1, len(train_data) + max_seq_id + 1))
# train_data["seq_id"] = train_index

# #rename pre fusion
# thermo_data.rename(columns = {'dtm':'tm', "gid":"group"}, inplace = True)

# # separate groups
# max_group_train = max(train_data["group"])
# thermo_data["group"] = thermo_data["group"] + max_group_train

# # filter invalid
# thermo_data = thermo_data[thermo_data["protein_sequence"] != "invalid"]
# thermo_data.dropna(axis=0, inplace=True)

# # fuse interesting variable
# train_tot = pd.concat([thermo_data[["protein_sequence", "group", "seq_id"]], train_data[["protein_sequence", "group", "seq_id"]]])
# train_tot.to_csv("train_tot.csv")
# save_data_as_fasta(test_data, "test.fasta", "seq_id")
# save_data_as_fasta(train_tot, "train_tot.fasta", "seq_id")



# Using train all

test_data = pd.read_csv("./test.csv")
save_data_as_fasta(test_data, "test.fasta", "seq_id")
train_data  = pd.read_csv("./train_all.csv", index_col=0)
save_data_as_fasta(train_data, "train.fasta")
