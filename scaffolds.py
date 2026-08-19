import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold

df = pd.read_csv("data/CHEMBL240_activities.csv")
smi = df["canonical_smiles"].iloc[0]

# -----------------------------------------------
# Initial test
# print(smi)

mol = Chem.MolFromSmiles(smi)
scaf = MurckoScaffold.GetScaffoldForMol(mol)
# print(Chem.MolToSmiles(scaf))
Draw.MolToFile(mol, "data/mol.png", size=(400, 400))

# -----------------------------------------------
# Actual function
uniq = df.drop_duplicates("canonical_smiles")
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
print(uniq["scaffold"].nunique(), "distinct scaffolds")
print("----------------")
print(
    uniq["scaffold"].value_counts().head(10).rename_axis("Top 10 Structures").toString()
)
