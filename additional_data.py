import json
import re
from io import StringIO
from operator import itemgetter

import pandas as pd
import requests
import sqlite3
from Bio import Seq, SeqIO

with open('databases/thermomutdb.json') as json_file:
    thermomut_db = json.load(json_file)

# filter database

# important columns
column_keep = ["dtm", "uniprot", "mutation_code", "effect"]
thermomut_db = [x for x in thermomut_db if all(
    v is not None for v in itemgetter(*column_keep)(x))]

# filter invalid uniport ids
thermomut_db = [x for x in thermomut_db if x["uniprot"].startswith("P")]

# filter by pH
thermomut_db_ph = [x for x in thermomut_db if x["ph"]]
thermomut_db_ph = [x for x in thermomut_db_ph if 5 < x["ph"] < 9]


# convert to dataframe
thermomut_df = pd.DataFrame(
    [itemgetter(*column_keep)(x) for x in thermomut_db], columns=column_keep)


# make groups
# get number of groups
thermomut_df.groupby("uniprot").ngroups
# make them and filter thoses with 4 or more proteins
thermomut_df = thermomut_df.groupby(
    "uniprot").filter(lambda x: len(x) >= 4)
# make a new column with the group id
thermomut_df['gid'] = (thermomut_df.groupby(
    ['uniprot']).cumcount() == 0).astype(int).cumsum()


# fetch sequences from uniprot
thermomut_df["protein_sequence"] = ""
base_url = "https://www.uniprot.org/uniprot/"

for index, row in thermomut_df.iterrows():
    currentUrl = base_url + thermomut_df.loc[index, "uniprot"] + ".fasta"
    response = requests.post(currentUrl)
    sequence = str(SeqIO.read(StringIO(''.join(response.text)), 'fasta').seq)
    mutation_code = thermomut_df.loc[index, "mutation_code"]
    match = re.findall(r"([a-z]+)([0-9]+)([a-z]+)", mutation_code, re.I)
    for m in match:
        pos = int(m[1])
        if pos > len(sequence):
            pos = len(sequence)
        if sequence[pos-1] == m[0]:
            sequence = sequence[:(pos-1)] + \
                m[2] + sequence[pos:]
            continue
        sequence = "invalid"
        break
    thermomut_df.loc[index, "protein_sequence"] = sequence

thermomut_df = thermomut_df[thermomut_df["protein_sequence"] != "invalid"]

# final dataframe
thermomut_df = thermomut_df[[
    "gid", "protein_sequence", "dtm", "effect"]].reset_index(drop=True)
