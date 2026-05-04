# Rainfall-induced microbial resuscitation reveals functional decoupling across biocrust succession

## Overview

This repository contains all code, workflows, and analysis scripts associated with:

> **CITATION**  
> *[Full paper title here]*  
> *Journal Name*, DOI: [XXXXXXX]

This project examines rainfall-induced microbial resuscitation dynamics in **cyanobacteria-dominated biological soil crusts (biocrusts)** across successional stages. Intact biocrust microcosms representing **early (L-BSC)** and late **(D-BSC)** stages were subjected to a 3 mm simulated precipitation event and incubated under light and dark conditions for 6 hr. Translationally active cells were identified using **BONCAT-FACS**, and both active and total (bulk DNA) communities were characterized via **16S rRNA gene amplicon sequencing**.

| Data type | Repository | Accession |
|---|---|---|
| Raw 16S gene sequences (.fastq.gz) | NCBI SRA | SRRXXXXXXX |
| BioProject | NCBI BioProject | PRJNA1458593 |
| BioSamples | NCBI BioSamples | SAMXXXXXXX |

## Code available
**Cutadap_resuscitation.sh** contains the pipeline for primer trimming using cutadapt.

**16S_resuscitation_pipeline.R** contains the 16S DADA2 Pipeline, 97% ASV clustering and taxonomy assignment

**16S_resuscitation_analysis.R** contains the 16S data analysis, including statistical analyses and plots.

**ANCOMBC-2_data_analysis.R** contains the ANCOMBC-2 output data analysis supporting figures 4 and 5.

## Citation

If you use this code or data, please cite:
