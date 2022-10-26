import numpy as np
import pandas as pd

# load datasets
train_df = pd.read_csv("train_grouped.csv", index_col=0)
thermomut_df = pd.read_csv("databases/thermomut_grouped.csv", index_col=0)

# prepare for merge
thermomut_df.drop(columns=["acc_id"], inplace=True)
thermomut_df.gid = thermomut_df.gid.add(train_df.gid.max())

# merge datasets
concat_train_df = pd.concat([train_df, thermomut_df], ignore_index=True)

# group by wildtype and redefine gid
concat_train_df['gid'] = concat_train_df.groupby("wildtype").ngroup().add(1)
concat_train_df.sort_values(by=["gid"], inplace=True)

concat_train_df.to_csv("train_all.csv")
