# Sorts through the dataframe to categorize molecules and sort them for usability
import numpy as np
import pandas as pd

df = pd.read_csv("data/CHEMBL240_activities.csv")
df_conf = pd.read_csv("data/CHEMBL240_confidence.csv")

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
pchembl = d[
    (d["parent_molecule_chembl_id"].notna())
    & (d["p_values"] > 3)
    & (d["potential_duplicate"] == 0)
]  # no missing ids


print("p-value summary")
print(
    pchembl.groupby("parent_molecule_chembl_id")["p_values"].agg(
        ["median", "std", "count"]
    )
)

pchembl = pchembl.drop_duplicates(
    subset=[
        "parent_molecule_chembl_id",
        "assay_chembl_id",
        "document_chembl_id",
        "standard_value",
    ]
)

pchembl_iterations = pchembl.groupby("parent_molecule_chembl_id")["p_values"].agg(
    ["median", "std", "count"]
)
pchembl_iterations = pchembl_iterations[pchembl_iterations["std"].notna()]
print(pchembl["p_values"].std(), "standard deviation of p-values")
print(pchembl_iterations)
print(pchembl_iterations["std"].median(), "median std of pchembl values")
print((pchembl_iterations["std"] == 0).sum(), "zeroes found in std of p-values' median")
print("-------------------------")

# print(
#     pchembl[
#         pchembl["parent_molecule_chembl_id"].isin(
#             pchembl_iterations[pchembl_iterations["std"] == 0].index
#         )
#     ][
#         [
#             "parent_molecule_chembl_id",
#             "assay_chembl_id",
#             "document_chembl_id",
#             "document_year",
#             "standard_value",
#             "potential_duplicate",
#         ]
#     ].sort_values(["parent_molecule_chembl_id", "document_chembl_id"])
# )

print("-------------------------")
# Confidence Score Analysis and Matching
print("Overall confidence score listing")
conf_groups = df_conf.groupby("confidence_score")["confidence_score"].count()

print(conf_groups.to_string(header=False))

print("---")
print("Confidence score matched to filtered data IDs")
conf_groups = df_conf[df_conf["assay_chembl_id"].isin(pchembl["assay_chembl_id"])]
# print(df["assay_chembl_id"].nunique())
print(
    conf_groups.groupby("confidence_score")["confidence_score"]
    .count()
    .to_string(header=False)
)
# print(conf_groups)
