# Extracts data for every entry of CHEMBL240 from the ChEMBL database and saves it to a CSV file
from pathlib import Path

import pandas as pd
from chembl_webresource_client.new_client import new_client

Path("data").mkdir(exist_ok=True)

COLS = [
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    #   "pchembl_value", # Not extracted until hypothesis
    "assay_chembl_id",
    "assay_type",
    "data_validity_comment",
    "potential_duplicate",
    "target_chembl_id",
    "document_chembl_id",
    "document_year",
]

acts = new_client.activity.filter(target_chembl_id="CHEMBL240", assay_type="B").only(
    COLS
)

df = pd.DataFrame(acts)
df.to_csv("data/CHEMBL240_activities.csv", index=False)
print(df.shape)
