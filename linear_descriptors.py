from pathlib import Path

import function
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

Path("data").mkdir(exist_ok=True)

# -----------------------------------------------
# Data import and filtering
# -----------------------------------------------
df = pd.read_csv("data/CHEMBL240_activities.csv")
df_conf = pd.read_csv("data/CHEMBL240_confidence.csv")
d = function.filtered_list(df, df_conf)

d = d.drop_duplicates(subset=["parent_molecule_chembl_id"])
d = d[d["canonical_smiles"].notnull()]

# -----------------------------------------------
# mol conversion with rdkit
# -----------------------------------------------
mol = d["canonical_smiles"].apply(Chem.MolFromSmiles)

# -----------------------------------------------
# Descriptor calculation
# -----------------------------------------------
MW = mol.apply(Descriptors.MolWt)
CLogP = mol.apply(Descriptors.MolLogP)
TPSA = mol.apply(Descriptors.TPSA)
HBD = mol.apply(Descriptors.NumHDonors)
Rotatable_Bonds = mol.apply(Descriptors.NumRotatableBonds)

# -----------------------------------------------
# CSV creation - 5 columns
# -----------------------------------------------
df = pd.DataFrame(
    {
        "parent_molecule_chembl_id": d["parent_molecule_chembl_id"],
        "MW": MW,
        "CLogP": CLogP,
        "TPSA": TPSA,
        "HBD": HBD,
        "Rotatable_Bonds": Rotatable_Bonds,
    }
)

df.to_csv("data/CHEMBL240_descriptors.csv", index=False)
