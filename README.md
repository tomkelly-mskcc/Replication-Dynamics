## Title
Modeling S. pombe Replication dynamics

## Project
The project implements a simple probabilistic model for analyzing DNA combing data and polymerase usage data.  for more information email tkelly@mskcc.org 

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
|CHROMOSOME_II_GENBANK	    |Chromosome-II-Anotation/ Schizosaccharomyces_pombe.ASM294v2.23.II.genbank|

## Parameter Scanning
The key replication parameters that define the dynamics of replication can be scanned over ranges of values by a set of recipes given in the file “src/DefaultParameterRecipes.java” 

## License
This project is licensed under the MIT License - see the LICENSE.md file for details

