#!/usr/bin/env python3
import blosum
import pandas as pd
import Levenshtein
import numpy as np
import matplotlib.pyplot as plt

DELETION_VALUE = -6
blo_mat = blosum.BLOSUM(90)

test_data = pd.read_csv("./test.csv")
bl90 = blosum.BLOSUM(90)
orig_prot = "VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK"
prot_pos = test_data[test_data["protein_sequence"] == orig_prot].index
test_data["new_order"] = range(0,len(test_data))
test_data.loc[prot_pos, "new_order"] = len(test_data)
test_data = test_data.sort_values("new_order").reset_index(drop="True").drop("new_order", axis=1)

edits = pd.Series(map(lambda x: Levenshtein.editops(orig_prot, x), test_data["protein_sequence"]))
len_edits = pd.Series(map(lambda x: len(x), edits))
edit_types = pd.Series(map(lambda x : x[0][0], edits[:len(edits) - 1]))
edit_types.value_counts()
set(len_edits) #only one at a time

def get_blosum_score_from_edit(lev_edit, mod_seq, orig_seq=orig_prot):
    if not lev_edit:
        return 0
    if lev_edit[0][0] == "delete":
        return DELETION_VALUE
    else:
        return blo_mat[orig_seq[lev_edit[0][1]] + mod_seq[lev_edit[0][2]]]

blo_scores = pd.Series(map(get_blosum_score_from_edit, edits, test_data.iloc[:, 1]))
test_data["scores"] = blo_scores

idiotic_submission = test_data[["seq_id", "scores"]].sort_values(by="scores", ascending=True)

order = list(range(idiotic_submission.shape[0]))
idiotic_submission["tm"] = order
idiotic_submission.drop("scores", axis=1, inplace=True)
idiotic_submission.to_csv("baseline.csv", index=False)

def get_blosum_score(mod_seq, orig_seq=orig_prot):
    lev_edit = Levenshtein.editops(orig_seq, mod_seq)
    if not lev_edit:
        return 0
    if lev_edit[0][0] == "delete":
        return DELETION_VALUE
    elif lev_edit[0][0] == "insert":
        return DELETION_VALUE
    else:
        res = 0
        for edit in lev_edit:
            res += blo_mat[orig_seq[edit[1]] + mod_seq[edit[2]]]
        return res

comp_matrix = np.empty((len(test_data), len(test_data)))
for i in range(len(comp_matrix)):
    for j in range(i+1, len(comp_matrix)):
        comp_matrix[i][j] = get_blosum_score(
            test_data.loc[i, "protein_sequence"],
            test_data.loc[j, "protein_sequence"]
        )


row_sums = np.sum(comp_matrix, axis=0)

less_idiotic_submission = pd.concat({"seq_id": test_data["seq_id"], "score": pd.Series(row_sums)}, axis=1)
less_idiotic_submission = less_idiotic_submission.sort_values(by="score", ascending=True)

plt.hist(less_idiotic_submission["score"])
plt.show()
order = list(range(less_idiotic_submission.shape[0]))
less_idiotic_submission["tm"] = order
less_idiotic_submission.drop("score", axis=1, inplace=True)
less_idiotic_submission.to_csv("baseline_2.csv", index=False)

