## Title
Modeling S. pombe Replication dynamics

## Project
The project implements a probabilistic model for analyzing DNA combing data and polymerase usage data.  for more information email tkelly@mskcc.org 

## Getting Started
The main method for analysis of combing data is in the file: “src/CombingExperimentSimulation.java”.

The main method for analysis of polymerase usage data is in the file:
 “src/ PuExperimentSimulation.java”

## Data files
The paths to required data files (and other global parameters) are defined as static variables in the file:
src/ ParameterSet.java

The value static variables holding the file paths in src/ ParameterSet.java should be changed to the locations of the following files:

|**variable name**        |**name of file or folder in this repository**|
|-------------             |---------------------------------------------|
|FA_SEQ_FILE              |Chromosome-II-Sequence/Chromosome-II_Sequence.fa|
|COMBING_DATA_FOLDER      |DNA-Combing-Data|
|BIN_COUNT_FILE           |PU-Seq-DataChromosome_II_bin_counts.csv|
|CHROMOSOME_II_GENBANK	    |Chromosome-II-Anotation/Schizosaccharomyces_pombe.ASM294v2.23.II.genbank|
|RIF1_INTERVAL_FILE       |Rif1-Sites/Rif1Intervals.txt|

The DNA cobing data is from  Kaykov A & Nurse P (2015) The spatial and temporal organization of origin firing during the S-phase of fission yeast. Genome Res 25(3):391-401.

The DNA polymerase usage data are from Daigaku Y, et al. (2015) A global profile of replicative polymerase usage. Nat Struct Mol Biol 22(3):192-198.

The DNA seqeunce and annotation of *S. pombe* chromosome II is from GENBANK.

## Parameter Scanning
The key replication parameters that define the dynamics of replication can be scanned over ranges of values by a set of recipes given in the file “src/DefaultParameterRecipes.java” 

## License
This project is licensed under the MIT License - see the LICENSE.md file for details

