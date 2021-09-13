# LBC_ARCHER

Repository for analysis of LBC ARCHER data.

# Longitudinal dynamics of clonal haematopoiesis identifies gene-specific fitness effects

Neil A. Robertson1+, Eric Latorre-Crespo1+, Maria Terradas-Terradas2,3, Alison C. Purcell2,3, Benjamin J Livesey1, Joseph A. Marsh1, Lee Murphy4, Angie Fawkes4, Louise MacGillivray4, Mhairi Copland2, Riccardo E. Marioni5, Sarah E. Harris6,7, Simon R. Cox6,7, Ian J. Deary6,7, Linus J. Schumacher8*, Kristina Kirschner2,3*, Tamir Chandra1*

1. MRC Human Genetics Unit, University of Edinburgh, Edinburgh, EH4 2XU, UK
2. Institute of Cancer Sciences, University of Glasgow, Glasgow, G61 1BD, UK
3. Cancer Research UK Beatson Institute, Glasgow, UK
4. Edinburgh Clinical Research Facility, University of Edinburgh, Edinburgh, EH4 2XU, UK
5. Centre for Genomic and Experimental Medicine, Institute of Genetics and Molecular Medicine, University of Edinburgh, Edinburgh, EH4 2XU, UK
6. Lothian Birth Cohorts, Department of Psychology, The University of Edinburgh, Edinburgh, UK
7. Department of Psychology, The University of Edinburgh, Edinburgh, UK
8. Centre for Regenerative Medicine, University of Edinburgh, Edinburgh, EH16 4UU, UK

\+ Equal Contribution

\* Correspondence to: linus.schumacher@ed.ac.uk, kristina.kirschner@glasgow.ac.uk, tamir.chandra@igmm.ed.ac.uk

![image](https://user-images.githubusercontent.com/4477113/123115304-eaa67780-d437-11eb-9adf-6a892a334b50.png)

The prevalence of clonal haematopoiesis of indeterminate potential (CHIP) in healthy individuals increases rapidly from age 60 onwards and has been associated with increased risk for malignancy, heart disease and ischemic stroke. CHIP is driven by somatic mutations in stem cells that are also drivers of myeloid malignancies. Since mutations in stem cells often drive leukaemia, we hypothesised that stem cell fitness substantially contributes to transformation from CHIP to leukaemia. Stem cell fitness is defined as the proliferative advantage over cells carrying no or only neutral mutations. It is currently unknown whether mutations in different CHIP genes lead to distinct fitness advantages that could form the basis for patient stratification. We set out to quantify the fitness effects of CHIP drivers over a 12 year timespan in older age, using longitudinal error-corrected sequencing data. We developed a new method based on drift-induced fluctuation (DIF) filtering to extract fitness effects from longitudinal data, and thus quantify the growth potential of variants within each individual. Our approach discriminates naturally drifting populations of cells and faster growing clones, while taking into account individual mutational context. We show that gene-specific fitness differences can outweigh inter-individual variation and therefore could form the basis for personalised clinical management.

Pre-print now online at: https://www.biorxiv.org/content/10.1101/2021.05.27.446006v3

# Repository overview

## General overview
This repository aims to develop tools to track and analyse the longitudinal
behaviour of genetic clones. The overarching goal is to infer the clonal
fitness or proliferative advantage of potentially malignant mutant clones.
The results of this analysis can be found in:
https://www.biorxiv.org/content/10.1101/2021.05.27.446006v3


This repository consists of several python modules associated with different
stages of the analysis process.

1. Transform data produced by ArcherDX into individual profiles with
trackable trajectories associated to genetic clones.

2. Create a filter to select clones presenting an abnormal behaviour.
    * Threshold filter: Selects variants based on a threshold of the clone size.
    * Neutral growth filter: Selects variants whose growth cannot be explained by
    stochastic drift of neutral clones.

3. Fit a model of exponential growth of stem cells to infer the fitness
advantage conferred by genetic variants to haematopoietic stem cells.

4. Quality control, statistical analysis and association with protein damage
prediction.

## Set up
### Data
This package has been designed to transform and analyse genetic data from the
Lothian Birth Cohort. To test this package ensure to transfer all supplementary
files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178936
to the Datasets folder.

### Python
A conda environment is included in this repository to facilitate the use of
this package.
'full_run.sh' script will install the conda environment and run all the
necessary python scripts to obtain all figures included in the article.

Alternatively, to install the environment run the following commands
(only valid for ubuntu):

    # extract and create environment from repository
    conda env create --file env/ARCH_env.yml --name ARCH_env
    # activate environment
    conda activate ARCH_env

    # add environment to jupyter notebooks
    conda install -c anaconda ipykernel
    python -m ipykernel install --user --name=ARCH_env

If the environment is already active, use 'analysis_run.sh' instead of
'full_run.sh'.

## Test
In order to facilitate testing, a less computationally heavy analysis is
implemented in 'full_run.sh' and 'analysis_run.sh' with a typical run time of
10 minutes.

To run the analysis as performed in the article, make the following
modifications:
- Modify line 132 of 'Scripts/full_NGF_run.py'
- Modify line 108 of 'Scripts/full_threshold_run.py'

to increase the number of random initialisations of the model fitting process.
In the article both iteration numbers were increased to 500, resulting in a
typical run time of 2h.

## Analysis exploration
To further explore the analysis of longitudinal trajectories, a series of
jupyter notebooks have been included in the 'Notebooks' folder.

# Repository maintenance
Code managed by:
+ https://github.com/elc08
+ https://github.com/neilrobertson
