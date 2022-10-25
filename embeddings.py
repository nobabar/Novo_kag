#!/usr/bin/env python3

import pandas as pd


# Save data as fasta ---------------------------
def save_data_as_fasta(data, output_file, col_id):
    with open(output_file, "w") as fileout:
        for _, prot in data.iterrows():
            fileout.write(f'>protein_id: {prot[col_id]}\n')
            fileout.write(f'{prot["protein_sequence"]}\n')


test_data = pd.read_csv("./test.csv")
train_data  = pd.read_csv("./train.csv")
thermo_data = pd.read_csv("./thermomut_grouped.csv", index_col=0)
max_seq_id = max(test_data["seq_id"])
thermo_index = list(range(max_seq_id + 1, len(thermo_data) + max_seq_id + 1))
thermo_data["seq_id"] = thermo_index

save_data_as_fasta(test_data, "test.fasta", "seq_id")
save_data_as_fasta(train_data, "train.fasta", "seq_id")
save_data_as_fasta(thermo_data, "thermo.fasta", "seq_id")
