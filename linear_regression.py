import function
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error


def lr_scaffold_split(Activities, Confidence, Descriptors, Runs):
    d = function.filtered_list(Activities, Confidence)

    d["scaffolds"] = d["canonical_smiles"].apply(function.scaffold)

    d = d.drop_duplicates(subset=["parent_molecule_chembl_id"])
    d = d.dropna(subset=["scaffolds"])

    total_count = len(d["parent_molecule_chembl_id"])
    scaffolds = d["scaffolds"].unique().tolist()

    d_group = d.groupby("scaffolds")
    d_dic = d_group["parent_molecule_chembl_id"].apply(list).to_dict()

    mae_list = []

    for i in range(Runs):
        print("Run", i)

        np.random.seed(i)  # Set seed for reproducibility
        np.random.shuffle(scaffolds)

        # -----------------------------------------------
        # train-test split
        # -----------------------------------------------
        train_ids = []

        for scaffold in scaffolds:
            x = d_dic[scaffold]

            train_ids.extend(x)
            if len(train_ids) / total_count >= 0.8:
                break

        # print(len(train_ids), "train molecules")
        # print(total_count, "total molecules")

        test_ids = d[~d["parent_molecule_chembl_id"].isin(train_ids)][
            "parent_molecule_chembl_id"
        ].tolist()

        # -----------------------------------------------
        # linear regression model
        # -----------------------------------------------
        # x axis
        Descriptors_train = Descriptors[
            Descriptors["parent_molecule_chembl_id"].isin(train_ids)
        ]

        train_x = Descriptors_train[["MW", "CLogP", "HBD", "Rotatable_Bonds", "TPSA"]]

        Descriptors_test = Descriptors[
            Descriptors["parent_molecule_chembl_id"].isin(test_ids)
        ]
        test_x = Descriptors_test[["MW", "CLogP", "HBD", "Rotatable_Bonds", "TPSA"]]

        # y axis
        p_values_train = Descriptors_train["p_values"]
        p_values_test = Descriptors_test["p_values"]

        model = LinearRegression()
        model.fit(train_x, p_values_train)
        predictions = model.predict(test_x)
        mae = mean_absolute_error(p_values_test, predictions)

        mae_list.append(mae)

    return mae_list


def lr_rando_split(Activities, Confidence, Descriptors, Runs):
    d = function.filtered_list(Activities, Confidence)

    d["scaffolds"] = d["canonical_smiles"].apply(function.scaffold)

    d = d.drop_duplicates(subset=["parent_molecule_chembl_id"])
    d = d.dropna(subset=["scaffolds"])

    mae_list = []

    id_list = d["parent_molecule_chembl_id"].tolist()

    for i in range(Runs):
        print("Run", i)

        np.random.seed(i)  # Set seed for reproducibility
        np.random.shuffle(id_list)

        train_ids = id_list[: round(len(id_list) * 0.8)]

        test_ids = id_list[round(len(id_list) * 0.8) :]

        # -----------------------------------------------
        # linear regression model
        # -----------------------------------------------
        # x axis
        Descriptors_train = Descriptors[
            Descriptors["parent_molecule_chembl_id"].isin(train_ids)
        ]
        train_x = Descriptors_train[["MW", "CLogP", "HBD", "Rotatable_Bonds", "TPSA"]]

        Descriptors_test = Descriptors[
            Descriptors["parent_molecule_chembl_id"].isin(test_ids)
        ]
        test_x = Descriptors_test[["MW", "CLogP", "HBD", "Rotatable_Bonds", "TPSA"]]

        # y axis
        p_values_train = Descriptors_train["p_values"]
        p_values_test = Descriptors_test["p_values"]

        model = LinearRegression()
        model.fit(train_x, p_values_train)
        predictions = model.predict(test_x)

        mae = mean_absolute_error(p_values_test, predictions)
        mae_list.append(mae)

    return mae_list
