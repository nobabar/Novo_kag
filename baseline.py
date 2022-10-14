#!/usr/bin/env python3
import blosum
import pandas as pd
import Levenshtein

DELETION_VALUE = -6
blo_mat = blosum.BLOSUM(90)

test_data = pd.read_csv("./test.csv")
bl90 = blosum.BLOSUM(90)
orig_prot = "VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK"
except_prot_pos = test_data[test_data["protein_sequence"] != orig_prot].index
prot_pos = test_data[test_data["protein_sequence"] == orig_prot].index

list_edit = []
len_edit = []
edits = pd.Series(map(lambda x: Levenshtein.editops(orig_prot, x), test_data["protein_sequence"]))
len_edits = pd.Series(map(lambda x: len(x), edits))
edit_types = pd.Series(map(lambda x : x[0][0], edits[except_prot_pos]))
edit_types.value_counts()
set(len_edits) #only one at a time

def get_blosum_score(lev_edit, mod_seq):
    if lev_edit[0][0] == "delete":
        return DELETION_VALUE
    else:
        return blo_mat[orig_prot[lev_edit[0][1]] + mod_seq[lev_edit[0][2]]]

blo_scores = pd.Series(map(get_blosum_score, edits[except_prot_pos], test_data.iloc[except_prot_pos, 1]))
blo_scores.loc[1168.5] = 0 #comment faire autrement aaargh
blo_scores = blo_scores.sort_index().reset_index(drop=True)
test_data["scores"] = blo_scores

idiotic_submission = test_data[["seq_id", "scores"]].sort_values(by="scores", ascending=False)

order = list(range(idiotic_submission.shape[0]))
idiotic_submission["tm"] = order
idiotic_submission.drop("scores", axis=1, inplace=True)
