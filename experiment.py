import math

import function
import linear_regression
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df_act = pd.read_csv("data/CHEMBL240_activities.csv")
df_conf = pd.read_csv("data/CHEMBL240_confidence.csv")
df_desc = pd.read_csv("data/CHEMBL240_descriptors.csv")

d = function.filtered_list(df_act, df_conf)
noise_floor = function.noise_floor(d)

runs = 50

print("----------------------------------")
print("Scaffold split")
print("----------------------------------")
scaffold_split = linear_regression.lr_scaffold_split(d, df_desc, runs)
print(scaffold_split)
print(".--.-.-.-.--..-.-")
print("MAE mean:", np.mean(scaffold_split))
print("Standard Deviation:", np.std(scaffold_split, ddof=1))

print("----------------------------------")
print("Random split")
print("----------------------------------")
random_split = linear_regression.lr_rando_split(d, df_desc, runs)
print(random_split)
print(".--.-.-.-.--..-.-")
print("MAE mean:", np.mean(random_split))
print("Standard Deviation:", np.std(random_split, ddof=1))

# --------------------------------------
# Pairing results
# --------------------------------------
paired_dif = np.array(scaffold_split) - np.array(random_split)
print("----------------------------------")
print("Paired differences")
print("----------------------------------")
print(paired_dif)
print(".--.-.-.-.--..-.-")
avg = np.mean(paired_dif)
print("Paired difference mean:", avg)
SD = np.std(paired_dif, ddof=1)
print("Paired difference standard deviation:", SD)

print(".--.-.-.-.--..-.-")
SE = SD / math.sqrt(runs)
print("Standard Error:", SE)

confidence_interval = [avg + 2 * SE, avg - 2 * SE]
print("Confidence interval:", confidence_interval)
print("Noise floor:", noise_floor)

print("////////////////////////------------------////////////////////////")
print("////////////////////////------------------////////////////////////")
if noise_floor < confidence_interval[0] or -noise_floor > confidence_interval[1]:
    print("The difference in split analysis is significant")
else:
    print("The difference in split analysis is not significant")
print("////////////////////////------------------////////////////////////")
print("////////////////////////------------------////////////////////////")

# --------------------------------------
# Pairing Difference Figure
# --------------------------------------
plt.hist(paired_dif)
plt.xlabel("Paired Differences")
plt.ylabel("Frequency")
plt.title("Distribution of Paired Differences")
plt.axvline(
    x=confidence_interval[0], color="r", linestyle="--", label="Confidence Interval"
)
plt.axvline(x=confidence_interval[1], color="r", linestyle="--", label="_nolegend")

plt.axvline(x=noise_floor, color="g", linestyle="-", label="Noise Floor")
plt.axvline(x=-noise_floor, color="g", linestyle="-", label="_nolegend")

plt.legend()
plt.xlim(-0.25, 0.25)
plt.savefig("figures/paired_differences.png")
plt.show()
