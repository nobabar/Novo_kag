import pandas as pd
from collections import Counter
from operator import itemgetter
from statistics import mode, StatisticsError

train_data = pd.read_csv("./train.csv", index_col="seq_id")

# apply update on train data
train_updates = pd.read_csv("train_updates_20220929.csv", index_col="seq_id")

all_features_nan = train_updates.isnull().all("columns")

drop_indices = train_updates[all_features_nan].index
train_data = train_data.drop(index=drop_indices)

swap_ph_tm_indices = train_updates[~all_features_nan].index
train_data.loc[swap_ph_tm_indices, ["pH", "tm"]
               ] = train_updates.loc[swap_ph_tm_indices, ["pH", "tm"]]

train_data.to_csv("./train_updated.csv")

# filter data
# by size
train_data = train_data[(train_data["protein_sequence"].str.len() < 1474)
                        & (train_data["protein_sequence"].str.len() > 89)]
# by pH
train_data = train_data[(train_data["pH"] > 5) & (train_data["pH"] < 9)]

train_data = train_data.sort_index().reset_index(drop=True)


# INSERTION DELETION THRESHOLD
D_THRESHOLD = 1
# MIN GROUP SIZE
MIN_GROUP_SIZE = 4


def max_item_count(seq):
    d = Counter(seq)
    return max(d.items(), key=itemgetter(1))


def mode_wildtype(proteins):
    wildtype = []
    try:
        for i in range(len(proteins.iloc[0])):
            wildtype.append(mode([p[i] for p in proteins]))
        return "".join(wildtype)
    except StatisticsError:
        return ""


def get_wildtype(proteins, is_retry=False):
    # try with mode solution
    if not is_retry:
        wildtype = mode_wildtype(proteins)
        if wildtype:
            return wildtype

    # at least 1/3rd length consecutive string must match. Find max counts of starts, middles, and ends
    # This technically isn"t a guaranteed or precise algorithm, but it is fast
    k = len(proteins.iloc[0]) // 3

    starts = [p[:k] for p in proteins]
    middles = [p[k:2*k] for p in proteins]
    ends = [p[-k:] for p in proteins]

    # get the most common substring, and the count of that substring
    start = max_item_count(starts)
    middle = max_item_count(middles)
    end = max_item_count(ends)

    # reduce the proteins to the ones that match the most common substring
    if (start[1] >= middle[1]) and (start[1] >= end[1]) and (start[1] >= MIN_GROUP_SIZE):
        proteins = pd.Series([p for p in proteins if p[:k] == start[0]])
        assert (start[1] == len(proteins))
    elif (middle[1] >= end[1]) and (middle[1] >= MIN_GROUP_SIZE):
        proteins = pd.Series([p for p in proteins if p[k:2*k] == middle[0]])
        assert (middle[1] == len(proteins))
    elif end[1] >= MIN_GROUP_SIZE:
        proteins = pd.Series([p for p in proteins if p[-k:] == end[0]])
        assert (end[1] == len(proteins))
    else:
        return ""

    # try again the mode solution with the reduced list
    return mode_wildtype(proteins)


train_data["gid"] = -1
train_data["wildtype"] = ""
grp = 0

train_data["length"] = train_data.protein_sequence.str.len()
vc = train_data.length.value_counts()

for i, count in enumerate(vc):
    if count < MIN_GROUP_SIZE:
        break
    length = vc.index[i]
    is_retry = False
    # subset train data with same protein length
    # doesn't take account of deletion for now for simplicity
    tmp = train_data.loc[(train_data.length == length)
                         & (train_data.gid == -1)]

    # It is possible that the same length protein string might have multiple wildtypes
    # keep searching until we've found all of them
    while len(tmp) >= MIN_GROUP_SIZE:
        # Ignore Levenstein distance, which is overkill
        # Directly attempt to find wildtype
        # Drop duplicates for wildtype guesstimation
        proteins = tmp.protein_sequence.drop_duplicates()

        # Create most likely wildtype
        wildtype = get_wildtype(proteins, is_retry=is_retry)
        if wildtype == "":
            break

        # subset train data with same protein length plus minus D_THRESHOLD
        tmp = train_data.loc[(train_data.length >= length - D_THRESHOLD) &
                             (train_data.length <= length + D_THRESHOLD) &
                             (train_data.gid == -1)]

        for idx in tmp.index:
            p = train_data.loc[idx, "protein_sequence"]

            k = length // 3
            # Use fast method to guess that it is only a few mutation away.
            if (wildtype[:k] == p[:k]) or (wildtype[k:k*2] == p[k:2*k]) or (wildtype[-k:] == p[-k:]):
                train_data.loc[idx, "gid"] = grp
                train_data.loc[idx, "wildtype"] = wildtype

        if len(train_data.loc[train_data.gid == grp]) >= MIN_GROUP_SIZE:
            grp += 1
            is_retry = False
        else:
            train_data.loc[train_data.gid == grp, "gid"] = -1
            train_data.loc[train_data.gid == grp, "wildtype"] = ""
            # to avoid an infinite loop, break out if we"ve already failed last time
            if is_retry:
                break
            is_retry = True

        # Get ready for next loop
        tmp = train_data.loc[(train_data.length == length)
                             & (train_data.gid == -1)]

# remove rows with no gid
train_data = train_data.loc[train_data.gid != -1]

# remove possible duplicates
train_data.drop_duplicates(subset=["protein_sequence"], inplace=True)

# make sure all groups still have at least MIN_GROUP_SIZE
train_data = train_data.groupby(
    "gid").filter(lambda x: len(x) >= MIN_GROUP_SIZE)

train_data.sort_values(by=["gid"], inplace=True)

# rank normalize the Tm
train_data["dtm"] = (10 + train_data.groupby("gid")["tm"].rank(method="dense")
                     ) / (20 + train_data.groupby("gid")["gid"].transform(len))

# add unique sequence id and rename group id
train_data["seqid"] = "S" + train_data.index.astype(str)
train_data["gid"] = "G" + train_data.gid.astype(str)

# final dataframe
train_data = train_data[["seqid", "gid", "protein_sequence",
                         "dtm", "wildtype"]].reset_index(drop=True)

train_data.to_csv("./train_grouped.csv")
