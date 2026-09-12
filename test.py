import linear_regression
import pandas as pd

df_act = pd.read_csv("data/CHEMBL240_activities.csv")
df_conf = pd.read_csv("data/CHEMBL240_confidence.csv")
df_desc = pd.read_csv("data/CHEMBL240_descriptors.csv")

print(linear_regression.linear_regression(df_act, df_conf, df_desc, 3))
