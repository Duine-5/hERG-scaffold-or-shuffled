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
# print(smi)

mol = Chem.MolFromSmiles(smi)
scaf = MurckoScaffold.GetScaffoldForMol(mol)
# print(Chem.MolToSmiles(scaf))
Draw.MolToFile(mol, "data/mol.png", size=(400, 400))

# -----------------------------------------------
# Actual function
uniq = d.drop_duplicates("canonical_smiles")
uniq = uniq.copy()


def scaffold(smi):
    if not isinstance(smi, str):  # check for non strings
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:  # Checks for non valid chemistry
        return None
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))


uniq["scaffold"] = uniq["canonical_smiles"].apply(scaffold)

print(uniq["scaffold"].isna().sum(), "failed to parse")
print(len(uniq), "distinct molecules")
print(uniq["scaffold"].nunique(), "distinct scaffolds")
print("----------------")
print(
    uniq["scaffold"]
    .value_counts()
    .head(10)
    .rename_axis("Top 10 Structures")
    .to_string()
)
print("----------------")

counts = uniq["scaffold"].value_counts()
print((counts == 1).sum(), "compounds with unique scaffolds")
print(
    round(counts[counts == 1].sum() / len(uniq) * (100), 2),
    "% of compounds have unique scaffolds",
)

Path("figures").mkdir(exist_ok=True)

counts.value_counts().sort_index().plot(
    kind="bar", logy=True
)  # vertical axis  on a log scale (logy)
plt.xlabel("Compounds per Scaffold")
plt.ylabel("log(Number of Scaffolds)")
plt.savefig("figures/scaffold_sizes.png", dpi=150, bbox_inches="tight")
