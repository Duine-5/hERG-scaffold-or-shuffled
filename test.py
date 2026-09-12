import linear_regression
import numpy as np
import pandas as pd

df_act = pd.read_csv("data/CHEMBL240_activities.csv")
df_conf = pd.read_csv("data/CHEMBL240_confidence.csv")
df_desc = pd.read_csv("data/CHEMBL240_descriptors.csv")

runs = 50

print("----------------------------------")
print("Scaffold split")
print("----------------------------------")
scaffold_split = linear_regression.lr_scaffold_split(df_act, df_conf, df_desc, runs)
print(scaffold_split)
print(".--.-.-.-.--..-.-")
print("MAE mean:", np.mean(scaffold_split))
print("Standard Deviation:", np.std(scaffold_split))

print("----------------------------------")
print("Random split")
print("----------------------------------")
random_split = linear_regression.lr_rando_split(df_act, df_conf, df_desc, runs)
print(random_split)
print(".--.-.-.-.--..-.-")
print("MAE mean:", np.mean(random_split))
print("Standard Deviation:", np.std(random_split))

# --------------------------------------
# Pairing results
# --------------------------------------
paired_dif = np.array(scaffold_split) - np.array(random_split)
print("----------------------------------")
print("Paired differences")
print("----------------------------------")
print(paired_dif)
print(".--.-.-.-.--..-.-")
print("Paired difference mean:", np.mean(paired_dif))
print("Paired difference mean:", np.std(paired_dif))
