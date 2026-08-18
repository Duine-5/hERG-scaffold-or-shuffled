import pandas as pd

df = pd.read_csv("data/CHEMBL240_activities.csv")

print(df.shape)
print(
    df["standard_type"].value_counts().head(10)
)  # how many times each value appears (Descending)
print(
    df["standard_relation"].value_counts()
)  # how many times each value appears (Descending)
print(
    df["molecule_chembl_id"].nunique()
)  # counts distinct values, counts how many molecules instead of values

print("----------------------")

d = df[df["standard_type"] == "IC50"]
print("IC50 only", len(d))

d = d[d["standard_relation"] == "="]
print("exact only", len(d))

d = d[d["standard_units"] == "nM"]
print("nM only", len(d))

print("unique molecules", d["molecule_chembl_id"].nunique())
