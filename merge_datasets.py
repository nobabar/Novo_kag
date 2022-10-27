import pandas as pd

# load datasets
train_df = pd.read_csv("train_grouped.csv", index_col=0)
thermomut_df = pd.read_csv("databases/thermomut_grouped.csv", index_col=0)
fireprot_df = pd.read_csv("databases/fireprot_grouped.csv", index_col=0)

# merge datasets
concat_train_df = pd.concat([train_df, thermomut_df, fireprot_df],
                            ignore_index=True)

# group by wildtype and redefine gid
concat_train_df['gid'] = "G" + \
    concat_train_df.groupby("wildtype").ngroup().astype(str).str.zfill(3)
concat_train_df.sort_values(by=["gid"], inplace=True)

# reset index and seqid
concat_train_df.reset_index(drop=True, inplace=True)
concat_train_df['seqid'] = "S" + concat_train_df.index.astype(str).str.zfill(5)

concat_train_df.to_csv("train_all.csv")

wildtypes = concat_train_df.groupby(
    "gid").wildtype.unique().apply(lambda x: x[0])

wildtypes.to_csv("wildtypes.csv")
