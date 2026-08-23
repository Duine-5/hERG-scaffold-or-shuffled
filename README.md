# Main question:
Will a random test split inflate apparent model performance relative to a scaffold based split for CHEMBL240? 

# Filter
- 32,640 bioactivity entries
- 16,198 entries which utilised IC50
- 10,228 entries which found IC50 within the tested concentrations
- 10,219 entries which report values in nM
- 8,394 unique compounds (Note: 1 compound could not be parsed)

## Findings
From 8,393 compounds, 4,266 scaffolds were found; 3059 compounds had a unique scaffold (36.45% of the total compounds), meaning that 1,207 scaffolds were repeated more than once across all compounds.

# Methods
## P Value
Concentrations above -3 (log) were discarded as it is the normal screening cutoff for inactive compounds. Results changed minimally (representing the removal of outlier testing):
- median 5.3 -> 5.301
- mean + 0.02
- std 1.008 -> 0.976

### Duplicate removal
Standard deviation across the same parent molecule was analysed, multiple duplicates skewed standard deviation towards 0. Duplicates were filtered through filtering out Potential Duplicates marked by Chembl and by filtering out data which contained same Parent Molecule ID, Assay ID, Document ID, and Standard Value.
> Note: Utilising the 4 filtering parameters allows some probable duplicate data in the system, but reduces the probability that genuine data is lost. There are still 40 entries (out of 8,304) which contain the exact same Standard Value, though they are not necessarily duplicates.
Standard deviation was found to be 0.265