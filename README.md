# Adzuki Bean Genomics Analysis Pipelines

Automated and reproducible bioinformatics workflows for whole-genome sequencing (WGS), population genomics, and evolutionary analyses on high-performance computing (HPC) environments.

## Overview

This repository contains a collection of automated analysis pipelines developed during large-scale genomic studies of adzuki bean (*Vigna angularis*) and related species.

The workflows are designed to enable reproducible, end-to-end analyses of next-generation sequencing (NGS) datasets on HPC systems with minimal manual intervention.

Most pipelines support one-command execution and integrate widely used bioinformatics tools with custom Perl, R, and shell scripts.

The repository was extensively used in population genomics studies, including:

> Chien, C.-C. et al. *Science* (2025). A single domestication origin of adzuki bean in Japan and the evolution of domestication genes.

---

## Key Features

* Automated end-to-end genomic analyses
* Reproducible HPC workflows
* Support for PBS, PBS Pro, and Slurm job schedulers
* One-command execution for complex analyses
* Modular pipeline structure
* Designed for large-scale WGS datasets
* Extensive use in published peer-reviewed studies

---

## Supported HPC Systems

* PBS
* PBS Pro
* Slurm

Most pipelines utilize custom job submission utilities:

```bash
qsub_subroutine.pl
create_job.pl
```

to automatically generate and submit jobs to HPC clusters.

---

## Pipeline Modules

| Directory | Description                                               |
| --------- | --------------------------------------------------------- |
| GATK4     | Variant calling workflow for whole-genome sequencing data |
| ADMIXTURE | Population structure analysis                             |
| GEMMA     | GWAS and mixed model analyses                             |
| MSMC      | Demographic history inference                             |
| Relate    | Genealogical and demographic analyses                     |
| TreeMix   | Population split and migration analyses                   |
| rehh      | Selection scan analyses (EHH/iHS/XP-EHH)                  |
| Beagle    | Genotype phasing and imputation                           |
| IQ-TREE   | Phylogenetic analyses                                     |
| STRUCTURE | Population assignment analyses                            |
| PCA       | Principal component analyses                              |
| FST       | Population differentiation analyses                       |

---

## Typical Workflow

```text
Raw FASTQ
   ↓
Quality Control
   ↓
Read Mapping
   ↓
Variant Calling (GATK4)
   ↓
Variant Filtering
   ↓
Population Genomics Analyses
   ├── PCA
   ├── ADMIXTURE
   ├── FST
   ├── TreeMix
   ├── MSMC
   ├── Relate
   ├── GWAS
   └── Selection Scan
   ↓
Visualization and Statistical Analyses
```

---

## Programming Languages

* Perl
* R
* Shell

---

## Dependencies

Examples include:

* GATK4
* BWA
* SAMtools
* BCFtools
* ADMIXTURE
* GEMMA
* MSMC2
* Relate
* TreeMix
* Beagle
* IQ-TREE
* PLINK
* VCFtools
* BEDTools
* R (various packages)

Please refer to each module directory for specific requirements.

---

## Citation

If you use these pipelines in your research, please cite:

Chien, C.-C. et al. (2025). *A single domestication origin of adzuki bean in Japan and the evolution of domestication genes*. Science.

---

## Author

Chih-Cheng Chien

National Agriculture and Food Research Organization (NARO), Japan

Research interests:

* Population genomics
* Evolutionary genomics
* Bioinformatics
* Genome editing
* Crop domestication

