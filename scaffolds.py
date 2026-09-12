from pathlib import Path

import function
import pandas as pd
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold

df = pd.read_csv("data/CHEMBL240_activities.csv")
df_conf = pd.read_csv("data/CHEMBL240_confidence.csv")
d = function.filtered_list(df, df_conf)

smi = d["canonical_smiles"].iloc[0]
# -----------------------------------------------
# Initial test
# -----------------------------------------------

# print(smi)

mol = Chem.MolFromSmiles(smi)
scaf = MurckoScaffold.GetScaffoldForMol(mol)
# print(Chem.MolToSmiles(scaf))
Draw.MolToFile(mol, "data/mol.png", size=(400, 400))

# -----------------------------------------------
# Scaffold generation
# -----------------------------------------------

d = d.drop_duplicates(subset=["parent_molecule_chembl_id"])

d["scaffold"] = d["canonical_smiles"].apply(function.scaffold)

# print(
#     d[d["scaffold"].isna()][
#         ["assay_chembl_id", "scaffold", "parent_molecule_chembl_id"]
#     ]
# )

print(d["scaffold"].isna().sum(), "failed to parse")
print(len(d), "distinct molecules")
print(d["scaffold"].nunique(), "distinct scaffolds")
print("----------------")
print(
    d["scaffold"].value_counts().head(10).rename_axis("Top 10 Structures").to_string()
)
print("----------------")

counts = d["scaffold"].value_counts()
print((counts == 1).sum(), "compounds with unique scaffolds")
print(
    round(counts[counts == 1].sum() / len(d[d["scaffold"].notnull()]) * (100), 2),
    "% of compounds have unique scaffolds",
)

Path("figures").mkdir(exist_ok=True)

counts.value_counts().sort_index().plot(
    kind="bar", logy=True
)  # vertical axis  on a log scale (logy)
plt.xlabel("Compounds per Scaffold")
plt.ylabel("log(Number of Scaffolds)")
plt.savefig("figures/scaffold_sizes.png", dpi=150, bbox_inches="tight")
