# Main question:
Will a random test split inflate apparent model performance relative to a scaffold based split for CHEMBL240? 

# Filter
- 32,640 bioactivity entries
- 16,198 entries which utilised IC50
- 10,228 entries which found IC50 within the tested concentrations
- 10,219 entries which report values in nM
- 8,394 -> 8,258 unique compounds (Note: 1 compound could not be parsed, change ocurred after data validity comments filter). 
- 8,202 unique molecules after parent molecule filter

## Scaffolds
From 8,232 compounds, 4,188 scaffolds were found; 3001 compounds had a unique scaffold (36.46% of the total compounds), meaning that 1,187 scaffolds were repeated more than once across all compounds.

# Methods
## P Value
Concentrations above -3 (log) were discarded after looking at the data, decision was formualted with basis on the normal screening cutoff for inactive compounds. Results changed minimally (representing the removal of outlier testing):
> Note: Difference between distinct compounds
- median 5.3 -> 5.301 -> 5.300 (after vailidity comments filter)
- mean + 0.02 above median -> 0.196 above median (after filters)
- std 1.008 -> 0.976 -> 0.950 (after removal of duplicates) -> 0.933 (after vailidity comments filter)

### Duplicate removal
Standard deviation across the same parent molecule was analysed, multiple duplicates skewed standard deviation towards 0. Duplicates were filtered through filtering out Potential Duplicates marked by Chembl and by filtering out data which contained same Parent Molecule ID, Assay ID, Document ID, and Standard Value.
> Note: Utilising the 4 filtering parameters allows some probable duplicate data in the system, but reduces the probability that genuine data is lost. There are still 40 entries (out of 8,304) which contain the exact same Standard Value, though they are not necessarily duplicates.

Standard deviation within same parent molecule (spread) was found to be 0.256

## Confidence Filter
A confidence filter with a threshold of >=8 was analysed, though no assay for this target scores below 8. After curation 19 assay ids were matched with a confidence level of 8, and 1,339 were matched with a confidence level of 9. No data would have been removed as a result of the filter as nothing scores below 8.
