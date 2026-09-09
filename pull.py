# Extracts data for every entry of CHEMBL240 from the ChEMBL database and saves it to a CSV file
from pathlib import Path

import function
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
conf = new_client.assay.filter(target_chembl_id="CHEMBL240", assay_type="B").only(
    "assay_chembl_id", "confidence_score"
)

if Path("data/CHEMBL240_activities.csv").exists() == False:
    df = pd.DataFrame(acts)
    df.to_csv("data/CHEMBL240_activities.csv", index=False)
    print(df.shape)
else:
    print("Activities skipped")

if Path("data/CHEMBL240_confidence.csv").exists() == False:
    df_c = pd.DataFrame(conf)
    print(df_c.shape)
    df_c.to_csv("data/CHEMBL240_confidence.csv", index=False)
else:
    print("Confidence skipped")

df_activities = pd.read_csv("data/CHEMBL240_activities.csv")
df_confidence = pd.read_csv("data/CHEMBL240_confidence.csv")
d = function.filtered_list(df_activities, df_confidence)
d = d["parent_molecule_chembl_id"].tolist()

InChIKeys = new_client.molecule.filter(molecule_chembl_id__in=d).only(
    "molecule_chembl_id",
    "molecule_structures",
    "canonical_smiles",
    "standard_inchi_key",
)
if Path("data/CHEMBL240_inchikey.csv").exists():
    print("InChIKeys skipped")
else:
    df = pd.DataFrame(InChIKeys)[["molecule_chembl_id", "molecule_structures"]].dropna()
    df["SMILES"] = df["molecule_structures"].map(lambda x: x["canonical_smiles"])
    df["InChIKey"] = df["molecule_structures"].map(lambda x: x["standard_inchi_key"])

    df.to_csv("data/CHEMBL240_inchikey.csv", index=False)
    print(len(InChIKeys), "InChIKeys")

# Test for diference in assay type = B and total entires in Chhembl
# print(len(new_client.activity.filter(target_chembl_id="CHEMBL240", assay_type="B")))
