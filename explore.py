# Sorts through the dataframe to categorize molecules and sort them for usability
import numpy as np
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
# Differentiate parent molecule from salt-forms
d = df[df["standard_type"] == "IC50"]
print("IC50 only", len(d))

d = d[d["standard_relation"] == "="]
print("exact only", len(d))

d = d[d["standard_units"] == "nM"]
print("nM only", len(d))

print("unique molecules", d["molecule_chembl_id"].nunique())

print("----------------------")
print("compound record counts")
by_mol = d.groupby("molecule_chembl_id").size().value_counts()
by_par = d.groupby("parent_molecule_chembl_id").size().value_counts()

Mol_comparison_tbl = pd.DataFrame({"molecule": by_mol, "parent": by_par}).sort_index()
Mol_comparison_tbl = pd.concat(
    [Mol_comparison_tbl, Mol_comparison_tbl.sum().to_frame("Total").T]
)
print(Mol_comparison_tbl)

print("----------------------")
d["p_values"] = 9 - np.log10(d.standard_value)  # Converts nm to -log(M)
pchembl = d[d["parent_molecule_chembl_id"].notna()]  # no missing ids

print("p-value summary")
print(pchembl["p_values"].describe())
