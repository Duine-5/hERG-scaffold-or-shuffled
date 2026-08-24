# Sorts through the dataframe to categorize molecules and sort them for usability
import function
import pandas as pd

df = pd.read_csv("data/CHEMBL240_activities.csv")
df_conf = pd.read_csv("data/CHEMBL240_confidence.csv")

# ----------------------------------------------------------------------------------------------------
# Preliminary dataframe analysis
# ----------------------------------------------------------------------------------------------------
print("-- Df rows and columns--")
print(df.shape)

print("-- Top 3 Standard type of assay --")
print(
    df["standard_type"].value_counts().head(3).to_string(header=False)
)  # how many times each value appears (Descending)

print("-- Standar Relations --")
print(
    df["standard_relation"].value_counts().to_string(header=False)
)  # how many times each value appears (Descending)

print("-- Unique Molecule IDs --")
print(
    df["molecule_chembl_id"].nunique(), "IDs"
)  # counts distinct values, counts how many molecules instead of values

# ----------------------------------------------------------------------------------------------------
# Filters
# ----------------------------------------------------------------------------------------------------
# print("----------------------")
# # Differentiate parent molecule from salt-forms
# d = df[df["standard_type"] == "IC50"]
# print("IC50 only", len(d))

# d = d[d["standard_relation"] == "="]
# print("exact only", len(d))

# d = d[d["standard_units"] == "nM"]
# print("nM only", len(d))

# d = d[
#     d["data_validity_comment"].isna()
# ]  # Filter based on analysis for data_validity comments, discarding all commented data

d = function.filtered_list(df, df_conf)

print("unique molecules", d["molecule_chembl_id"].nunique())

# ----------------------------------------------------------------------------------------------------
# Comparison between parent molecule id and moleculte id for salt-form analysis
# Commented after analysis
# ----------------------------------------------------------------------------------------------------
# print("----------------------")
# print("compound record counts")
by_mol = d.groupby("molecule_chembl_id").size().value_counts()
by_par = d.groupby("parent_molecule_chembl_id").size().value_counts()

Mol_comparison_tbl = pd.DataFrame({"molecule": by_mol, "parent": by_par}).sort_index()
Mol_comparison_tbl = pd.concat(
    [Mol_comparison_tbl, Mol_comparison_tbl.sum().to_frame("Total").T]
)
# print(Mol_comparison_tbl)

# ----------------------------------------------------------------------------------------------------
# p-values analysis, filter | Spread analysis | Duplicate filtering
# ----------------------------------------------------------------------------------------------------
print("----------------------")
# d["p_values"] = 9 - np.log10(d.standard_value)  # Converts nm to -log(M)
# pchembl = d[
#     (d["parent_molecule_chembl_id"].notna())
#     & (d["p_values"] > 3)
#     & (d["potential_duplicate"] == 0)

print("p-value summary")
# Table for p-value count analysis
# print(
#     pchembl.groupby("parent_molecule_chembl_id")["p_values"].agg(
#         ["median", "std", "count"]
#     )
# )

# pchembl = pchembl.drop_duplicates(
#     subset=[
#         "parent_molecule_chembl_id",
#         "assay_chembl_id",
#         "document_chembl_id",
#         "standard_value",
#     ]
# )

pchembl_iterations = d.groupby("parent_molecule_chembl_id")["p_values"].agg(
    ["median", "std", "count"]
)
pchembl_iterations = pchembl_iterations[pchembl_iterations["std"].notna()]
print(f"std: {d['p_values'].std():.3f} of p-values")
print(f"median: {d['p_values'].median():.3f} of p-values")
print(f"mean: {d['p_values'].mean():.3f} of p-values")
print(".......")
# print(pchembl_iterations) # Table for count analysis and comparison with non-filtered version
print(f"median std: {pchembl_iterations['std'].median():.3f} of pchembl values")
print((pchembl_iterations["std"] == 0).sum(), "zeroes found in std of p-values' median")

# ----------------------------------------------------------------------------------------------------
# Table to analyse potential duplicates
# Verification completed, currently not needed
# ----------------------------------------------------------------------------------------------------
# print("-------------------------")
# print(
#     d[
#         d["parent_molecule_chembl_id"].isin(
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

# ----------------------------------------------------------------------------------------------------
# Confidence Score Analysis and Matching
# Verification completed, current code is not needed
# ----------------------------------------------------------------------------------------------------
# print("-------------------------")
# print("Overall confidence score listing")
# conf_groups = df_conf.groupby("confidence_score")["confidence_score"].count()

# print(conf_groups.to_string(header=False))

# print("---")
# print("Confidence score matched to filtered data IDs")
# conf_groups = df_conf[df_conf["assay_chembl_id"].isin(pchembl["assay_chembl_id"])]
# # print(df["assay_chembl_id"].nunique())
# print(
#     conf_groups.groupby("confidence_score")["confidence_score"]
#     .count()
#     .to_string(header=False)
# )
# print(conf_groups)

# ----------------------------------------------------------------------------------------------------
# Data Validity Analysis
# Filter decided based on data explored here
# ----------------------------------------------------------------------------------------------------
# print("-------------------------")
# print("Data validity Analysis")
# data_validity_grouped = d.groupby("data_validity_comment")[
#     "data_validity_comment"
# ].count()
# print(data_validity_grouped.to_string())

# data_validity_check = d[d["data_validity_comment"].notna()]
# print(data_validity_check)
# print(data_validity_check["p_values"].agg(["min", "max", "median"]).to_string())
print("-------------------------")
# print(d["parent_molecule_chembl_id"].nunique())
# print(d["parent_molecule_chembl_id"].sample().item())
print(d.shape)
print(d["parent_molecule_chembl_id"].nunique())
