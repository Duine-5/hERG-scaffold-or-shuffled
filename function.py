import numpy as np


def filtered_list(df_activities, df_confidence):
    # df = pd.read_csv("data/CHEMBL240_activities.csv")
    # df_conf = pd.read_csv("data/CHEMBL240_confidence.csv")

    d = df_activities[
        (df_activities["standard_type"] == "IC50")
        & (df_activities["standard_relation"] == "=")
        & (df_activities["standard_units"] == "nM")
        & (df_activities["data_validity_comment"].isna())
        & (df_activities["parent_molecule_chembl_id"].notna())
        & (df_activities["potential_duplicate"] == 0)
    ]

    d["p_values"] = 9 - np.log10(d.standard_value)
    d = d[(d["p_values"] > 3)]

    d = d.drop_duplicates(
        subset=[
            "parent_molecule_chembl_id",
            "assay_chembl_id",
            "document_chembl_id",
            "standard_value",
        ]
    )

    d = d[
        d["assay_chembl_id"].isin(
            df_confidence[df_confidence["confidence_score"] >= 8]["assay_chembl_id"]
        )
    ]

    repetitions = d.groupby("parent_molecule_chembl_id").size()
    group_difference = d[
        d["parent_molecule_chembl_id"].isin(repetitions[repetitions >= 2].index)
    ].groupby("parent_molecule_chembl_id")
    min_max_p = group_difference["p_values"].agg([min, max])
    dif = min_max_p["max"] - min_max_p["min"]
    d = d[~d["parent_molecule_chembl_id"].isin(dif[dif > 1.36].index)]

    return d
