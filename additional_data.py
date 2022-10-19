import json
import re
from io import StringIO
from operator import itemgetter

import pandas as pd
import requests
from Bio import Entrez, SeqIO

with open('databases/thermomutdb.json') as json_file:
    thermomut_db = json.load(json_file)

# filter database

column_keep = ["dtm", "uniprot", "PDB_wild", "mut_count", "effect"]
thermomut_db = [x for x in thermomut_db if all(
    v is not None for v in itemgetter(*column_keep)(x))]


# entries with missing values in the following fields are removed
thermomut_db = [x for x in thermomut_db if all(
    itemgetter("ph", "length")(x))]
# by size
thermomut_db = [x for x in thermomut_db if 100 < x["length"] < 1000]
# by pH
thermomut_db = [x for x in thermomut_db if 5 < x["ph"] < 9]


thermomut_df = pd.DataFrame(
    [itemgetter(*column_keep)(x) for x in thermomut_db], columns=column_keep)

thermomut_df.groupby("uniprot").ngroups
thermomut_df.groupby("PDB_wild").ngroups

thermomut_df = thermomut_df.groupby(
    "uniprot").filter(lambda x: len(x) >= 4)

thermomut_df['gid'] = (thermomut_df.groupby(
    ['uniprot']).cumcount() == 0).astype(int).cumsum()

# fetch sequences from uniprot
thermomut_df["protein_sequence"] = ""
base_url = "https://www.uniprot.org/uniprot/"

for index, row in thermomut_df.iterrows():
    currentUrl = base_url + thermomut_df.loc[index, "uniprot"] + ".fasta"
    response = requests.post(currentUrl)
    Seq = list(SeqIO.parse(StringIO(''.join(response.text)), 'fasta'))[0]

    mutation_code = thermomut_df.loc[index, "mutation_code"]
    match = re.findall(r"([a-z]+)([0-9]+)([a-z]+)", mutation_code, re.I)

    for m in match:
        if Seq.seq[int(m[1])] == m[0]:
            Seq.seq = Seq.seq[:int(m[1])-1] + \
                m[2] + Seq.seq[int(m[1]):]

    thermomut_df.loc[index, "protein_sequence"] = str(Seq.seq)

# final dataframe
thermomut_df = thermomut_df[[
    "gid", "protein_sequence", "dtm", "effect"]].reset_index(drop=True)
