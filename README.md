# Experimental design
## Hypothesis
Random training data sets inflates linear regression model performance for pChembl relative to a scaffold-based split data set. 
> Note: drafted during commit f168b82 2026-08-27, edits made after incur corrections
## Tested Descriptors
Target will be CHEMBL240, with data extracted from CHEMBL 37
5 descriptors will be:
- MW ([Molecular Weight] Greater contact surface typically allows for greater potency)
- CLogP ([Calculated logP] Binding tends to work through hydrophobic interactions)
- TPSA ([Topological Polar Surface Area] A greater TPSA is an indicator of solubility and hydrophylic interactions, which may work against binding)
- HBD ([Hydrogen Bond Donor count] More N-H and O-H groups allows for binding)
- Rotatable Bonds (Entropy creates friction for high flexibility molecules to bind)

> Note: a 5 descriptor linear regression is expected to be a weakindicator for trends across drug-like chemistry, for that reason, this experiment relates specifically to a comparison in behaviour rather than a direct affinity score.
## Scaffold Framework
The core molecular backbone will be determined via Bemis-Murcko framework (Murcko Scaffold).
## Test set and Repeats
> Data will utilise random seeds for both distribution and repeats
The data will be split into 80% training data, 20% testing data. 
50 repeats will be performed 
## Success metrics
Gaps below the noise floor of 0.256 log units will be noted as uninterpretable. (See Duplicate Removal for explanation)
To accurately assess the model's performance, MAE was chosen as the success evaluation metric. 
## Statistical analysis
An alpha of 0.05, corresponding to a 95% will be utilised as the width of the paired interval.
A MAE paired interval between the random and scaffold split will determine the model's accuracy.

# Filter
- 41,078 bioactivity entries
- 32,640 binding assay entries
- 16,198 entries which utilised IC50
- 10,228 entries which found IC50 within the tested concentrations
- 10,219 entries which report values in nM
- 8,394 -> 8,258 unique compounds 
- 8,202 unique parents
- 8,200 molecule structures (CHEMBL5869391 and CHEMBL5784130 had no structure)

## Scaffolds
From 8,202 compounds, 4,187 scaffolds were found; 3014 compounds had a unique scaffold (36.76% of the total compounds), meaning that 1,173 scaffolds were repeated more than once across all compounds.

Compounds with no scaffold were dropped from the analysis as they cannot be parsed and do not represent a significant sample.

# Methods
## P Value
> Note: Outcome values were computed and inspected for curation, but were never measured against any predictor
Concentrations above -3 (log) were discarded after looking at the data, decision was formualted with basis on the normal screening cutoff for inactive compounds. Results changed minimally (representing the removal of outlier testing):
> Note: Difference between distinct compounds
- median 5.3 -> 5.301 -> 5.300 (after vailidity comments filter)
- mean + 0.02 above median -> 0.196 above median (after filters)
- std 1.008 -> 0.976 -> 0.950 (after removal of duplicates) -> 0.933 (after vailidity comments filter)

### Duplicate removal
Standard deviation across the same parent molecule was analysed, multiple duplicates skewed standard deviation towards 0. Duplicates were filtered through filtering out Potential Duplicates marked by Chembl and by filtering out data which contained same Parent Molecule ID, Assay ID, Document ID, and Standard Value.
> Note: Utilising the 4 filtering parameters allows some probable duplicate data in the system, but reduces the probability that genuine data is lost. There are still 40 entries (out of 8,304) which contain the exact same Standard Value, though they are not necessarily duplicates.

Standard deviation within same parent molecule (spread) was found to be 0.256, which was determined as the noise floor.

## Confidence Filter
A confidence filter with a threshold of >=8 was analysed, though no assay for this target scores below 8. After curation 19 assay ids were matched with a confidence level of 8, and 1,339 were matched with a confidence level of 9. No data would have been removed as a result of the filter as nothing scores below 8.
