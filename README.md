# GCS
The genomic complexity score (GCS) is able to enumerate the whole-genome complexity based on its respective copy number profile using SNP array or WGS data.

## Rational

In a diploid genome (2N) the log2 ratio is typically equal to 0, and it is graphically represented by probes or bins that cluster around 0 (Figure 1A). A deviation from this situation means that there is a loss or a gain of a chromosome or pieces of it, resulting in an increase of the complexity of the genome (Figure 1B). Assuming that there is a correlation between a low genome complexity and few copy number alterations genome-wide, and vice versa for a highly complex genome, then the more complex a genome is, the more the probes or the bins should deviate from a diploid copy number profile state.

![Fig_GCS](https://github.com/user-attachments/assets/b597b76a-ea2a-4dfe-b325-570c729da0a0)

Figure 1. Genomic complexity score (GCS) calculation.
Examples of complexity for chromosome 17. The blue rectangle marks the range where the probes cluster and the dashed lines are the boundaries based on the two values obtained at +1.96 and -1.96 standard deviations. (A) The graph shows a diploid situation with a low complex genome profile where the probes cluster close to 0. (B) The p arm of chromosome 17 is subjected to big changes of log2 intensity representing amplifications and copy number changes.

Taking all of these aspects into consideration, if we suppose that the mean of the log2 values of bins or probes of a chromosome without copy number changes is equal to 0, and same for the standard deviation, then there will also be a spread of the log2 when we have a change in the copy number profile. To assess how complex a genome is, based on the copy number profile, the mean and standard deviation of the log2 value of all genomic bins or individual probes were calculated for each chromosome.

The GCS can be correlated to e.g., the age of the patient (https://doi.org/10.1016/j.labinv.2023.100283 and https://doi.org/10.1002/path.6219) or to the number of chromosomes affected by a copy number changes, depending on the research question at hand. It can also be used to see differences in cases affected by a specific mutation, e.g. _TP53_.

## Getting Started

GCS is implemented in R and python. Before you can run this repository, make sure to install and set up:

- Git. `https://github.com/git-guides/install-git`
- Conda. `Install it through miniconda or anaconda  `
- R. 

## Project setup

### Installation

Clone the GCS repository fro the Github Page and execute the following to generate a conda environment with all you need to run the GCS

```
git clone https://github.com//ValeriaDifilippo/GCS
cd GCS
conda env create -f bin/config/environment.yml
conda activate GCS
```
or dowload manually the repository from the Git up page and then

```
cd GCS
conda env create -f bin/config/environment.yml
conda activate GCS
```

Create inputs and results directories in the NAFuse main folder 
```
mkdir -p input_files 
mkdir -p results
```

**Project folder structure**

The project is organized into the following main directories:

- _bin/_
This directory contains all the source codes and scripts used to run GCS and the jupyter notebook. In the _config/_ subdirectory you can find the file to set up your environment.

- _input_files/_
Contains the input files in `.txt` format for the analysis. 

- _results/_
This directory contains subdirectorie for each case with results generated by running GCS.

## Usage

1. Make sure to extract the columns needed for the analysis and to rename them like the example below. The input files need to be in a  tab delimited format and saved as `.txt` in the _input_files/_ directory.
   
**Columns** | **Explanation**
-- | --
Chromosome | 1 to 22 (autosome)
Position | The position of the probes or bins
LogR | Depending on the sequencing tecnique and the algorithm that you are using the column might be named differently 

2. Run The main core

```
  Rscript bin/GCS.R /path/to/GCS/
```


The final outputs of GCS will be saved in the _results/_ directory creating a subdirectory for each case. The tool will produce 

- `sample1_GCS.txt`
The **final genomic complecity score** given by the mean of the GCSs in the `sample1_all_chromosomes_results.txt`. (Recommanded to use this)

- `sample1_all_chromosomes_results.txt`
A file containing the information for each chromosome.

**Columns** | **Explanation**
-- | --
Case | Sample ID 
GCS | The genomic complecity score that represent the complexity of the genome for each chromosome.
Chr | 1 to 22 (autosome)
SD_pos | Value of the postive standard deviation
SD_neg | Value of the negative standard deviation

- `sample1_whole_genome_results.txt`
A file containing the genomic complecity score of the whole genome.

- _`plots/`_
A subdirectory containing the plot for each chromosome.
  
<img src="https://github.com/user-attachments/assets/f3d947d8-8498-4cdb-aa42-833c344b287c" alt="Untitled-1" width="50%">

## Additional feature

In the _bin/_ directory there is the jupyter notebook `GCS_plot.ipynb`. The code can plot the cases previosly analyzed in a dynamically way. 

```
jupyter-notebook bin/GCS_plot.ipynb
```

## Citations

Please if you use the GCS, cite:

Valeria Difilippo, Karim H. Saba, Emelie Styring, Linda Magnusson, Jenny Nilsson, Michaela Nathrath, Daniel Baumhoer, Karolin H. Nord, Osteosarcomas With Few Chromosomal Alterations or Adult Onset Are Genetically Heterogeneous, Laboratory Investigation, Volume 104, Issue 1, 2024, 100283, ISSN 0023-6837, https://doi.org/10.1016/j.labinv.2023.100283.
