# Longitudinal analysis of variants associated with ARCH

# General overview
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
Lothian Birth Cohort. To test this package ensure to transfer all files from
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178936
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
