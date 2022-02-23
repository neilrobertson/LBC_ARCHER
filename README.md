# LBC_ARCHER

Repository for analysis of LBC ARCHER data.

# Longitudinal dynamics of clonal haematopoiesis identifies gene-specific fitness effects

Neil A. Robertson1+, Eric Latorre-Crespo1+, Maria Terradas-Terradas2,3, Jorge Lemos-Portela4 Alison C. Purcell2,3, Benjamin J. Livesey1, Robert F. Hillary5, Joseph A. Marsh1, Lee Murphy6, Angie Fawkes6, Louise MacGillivray6, Mhairi Copland7, Riccardo E. Marioni5, Sarah E. Harris8,9, Simon R. Cox8,9, Ian J. Deary8,9, Linus J. Schumacher4#+, Kristina Kirschner2,3#+, Tamir Chandra1#+
 
1. MRC Human Genetics Unit, University of Edinburgh, Edinburgh, EH4 2XU, UK
2. Institute of Cancer Sciences, University of Glasgow, Glasgow, G61 1BD, UK
3. Cancer Research UK Beatson Institute, Glasgow, UK
4. Centre for Regenerative Medicine, University of Edinburgh, Edinburgh, EH16 4UU, UK
5. Centre for Genomic and Experimental Medicine, Institute of Genetics and Cancer, University 
6. Edinburgh Clinical Research Facility, University of Edinburgh, Edinburgh, EH4 2XU, UK
7. Paul O’Gorman Leukaemia Research Centre, Institute of Cancer Sciences, University of Glasgow, Gartnavel General Hospital, 1053 Great Western Road, Glasgow, G12 0YN, UK
of Edinburgh, Edinburgh, EH4 2XU, UK
8. Lothian Birth Cohorts, Department of Psychology, The University of Edinburgh, Edinburgh, UK
9. Department of Psychology, The University of Edinburgh, Edinburgh, UK

\+ Equal Contribution

\* Correspondence to: linus.schumacher@ed.ac.uk, kristina.kirschner@glasgow.ac.uk, tamir.chandra@igmm.ed.ac.uk

<img width="807" alt="image" src="https://user-images.githubusercontent.com/4477113/155323166-42d30bfe-d1dd-47b8-aa91-f404fcdc9f17.png">

The prevalence of clonal haematopoiesis of indeterminate potential (CHIP) in healthy individuals increases rapidly from age 60 onwards and has been associated with increased risk for malignancy, heart disease and ischemic stroke. CHIP is driven by somatic mutations in stem cells that are also drivers of myeloid malignancies. Since mutations in stem cells often drive leukaemia, we hypothesised that stem cell fitness substantially contributes to transformation from CHIP to leukaemia. Stem cell fitness is defined as the proliferative advantage over cells carrying no or only neutral mutations. It is currently unknown whether mutations in different CHIP genes lead to distinct fitness advantages that could form the basis for patient stratification. We set out to quantify the fitness effects of CHIP drivers over a 12-year timespan in older age, using longitudinal error-corrected sequencing data. Two key results support the possibility for individualised clinical monitoring of CHIP: (i) We developed a new filtering method to extract fitness effects from longitudinal data using Bayesian inference, while taking into account individual mutational context and co-occurrence of mutations, and thus quantify the growth potential of variants within each individual. (ii) We show that gene-specific fitness differences can outweigh inter-individual variation and therefore could form the basis for personalised clinical management in the future.

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
    * Likelihood-based Filter for Time-series data (LIFT): Selects variants using 
    Bayesian model comparison.

3. Use Bayesian model comparison to select the most likely clonal structure 
using a model of exponential growth of stem cells and infer the fitness
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

## Test
In order to facilitate testing, a less computationally heavy analysis  can be
implemented modifying the resolution of the implementation of the Bayesian model
comparison.

Set global parameters:
- Modify lines 40-45 in 'Scripts/LiFT_clonal_fit.py'
- Modify lines 39- 44 in 'Scripts/threshold_clonal_fit.py'

Currently, the full implementation of the code runs in ~12 hours on 60 cores.

## Analysis exploration
To further explore the analysis of longitudinal trajectories, a series of
jupyter notebooks have been included in the 'Notebooks' folder.

# Repository maintenance
Code managed by:
+ https://github.com/elc08
+ https://github.com/neilrobertson
