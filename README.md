# ccle_processing
What you need to process the Quarterly DepMap-Omics releases from Terra

## Instalation

/!\ this repository needs some important data and code from the [JKBio Library](https://www.github.com/jkobject/JKBio)
/!\ you would need the approriate R packages and python packages:
- for R packages, a loading function contains all required ones (in src's R file)
- for Python *no python requirements created yet*

## File structure

there is for now 3 computation pipeline for depmap omics:
- expression
- mutations
- copy number

each:
gets data from buckets, 
updates the TSVs on Terra, 
compute the results for each, 
QC them and do some more filtering,
Uploads them to taiga.

data/ contains important information used for processing
src contains the location of function files
