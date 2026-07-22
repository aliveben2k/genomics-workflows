# Genomics Analysis Pipelines

Automated and reproducible bioinformatics workflows for whole-genome sequencing (WGS), population genomics, and evolutionary genomics analyses on High-Performance Computing (HPC) environments.
Some pipelines have Nextflow or Nextflow+Python version, please check the folder "Nextflow_version" first.

---

## Overview

This repository contains a collection of automated bioinformatics workflows developed during large-scale genomic studies of adzuki bean (*Vigna angularis*) and related species.

The core component of this repository is an automated GATK-based variant calling framework designed for large-scale next-generation sequencing (NGS) analyses on HPC systems. The pipeline enables end-to-end processing of raw sequencing data, allowing users to generate high-quality VCF files directly from raw FASTQ files using a single command line.

The automated GATK workflow supports:

* Whole-genome sequencing (WGS)
* Whole-exome sequencing (WES)
* RAD-seq / DArT-seq data
* GPU-accelerated variant calling using NVIDIA Parabricks
* Incremental sample integration
* Joint genotyping of hundreds to thousands of samples
* Automated HPC job submission and dependency management
* Germline CNV analysis

In addition to variant calling, this repository also provides a collection of downstream population genomics and evolutionary genomics workflows for population structure analysis, demographic inference, phylogenetic reconstruction, selection scans, genotype imputation, and genome-wide association studies (GWAS).

Most workflows support one-command execution and integrate widely used bioinformatics tools with custom Perl, R, and shell scripts to enable reproducible, end-to-end genomic analyses on HPC systems with minimal manual intervention.

The repository was extensively used in population genomics studies, including:

> Chien, C.-C., Seiko, T., Muto, C., Ariga, H., Wang, Y.-C., Chang, C.-H., Sakai, H., Naito, K., & Lee, C.-R. (2025). *A single domestication origin of adzuki bean in Japan and the evolution of domestication genes*. **Science, 388**(6750), eads2871.

---

## Key Features

* Automated end-to-end genomic analyses
* Reproducible HPC workflows
* One-command execution from raw FASTQ files to final VCF files
* Support for PBS, PBS Pro, and Slurm job schedulers
* Configuration-driven deployment across multiple HPC environments
* Automated HPC job submission and dependency management
* Modular pipeline architecture
* GPU-accelerated variant calling support
* Designed for large-scale WGS datasets
* Extensively used in peer-reviewed publications

---

## Supported HPC Systems

* PBS
* PBS Pro
* Slurm
* SGE/UGE

Most pipelines utilize custom job submission utilities:

```bash
qsub_subroutine.pl
create_job.pl
```

to automatically generate and submit jobs to HPC clusters.

The workflow framework supports configuration-driven deployment across multiple HPC environments through a user-defined configuration file (`qsub_server.conf`).

The server configuration file can be provided in one of the following ways (listed in order of priority):

### 1. Explicitly specify the configuration file

```bash
perl create_job_v2.pl -cj_server /path/to/qsub_server.conf
```

### 2. Place the configuration file in the user's software directory

```text
$HOME/software/qsub_server.conf
```

### 3. Place the configuration file in the current working directory

```text
./qsub_server.conf
```

This configuration mechanism allows the same workflow framework to be easily deployed across different HPC clusters without modifying the source code.

---

## Pipeline Modules

### Variant Calling

| Directory              | Description                                                                                                                                            |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| GATK4                  | Automated Perl-based variant calling workflow for whole-genome sequencing data                                                                         |
| nf-core-gatkpopulation | Nextflow implementation supporting CPU/GPU germline calling, WGS/WES/RAD/amplicon data, population joint genotyping, and multiple HPC schedulers |

### Downstream Population Genomics and Evolutionary Genomics Analyses

| Directory | Description                                                     |
| --------- | --------------------------------------------------------------- |
| ADMIXTURE | Population structure analysis                                   |
| GEMMA     | Genome-wide association studies (GWAS) and mixed model analyses |
| MSMC      | Demographic history inference                                   |
| Relate    | Genealogical and demographic analyses                           |
| TreeMix   | Population split and migration analyses                         |
| rehh      | Selection scan analyses (EHH/iHS/XP-EHH)                        |
| Beagle    | Genotype phasing and imputation                                 |
| IQ-TREE   | Phylogenetic analyses                                           |
| STRUCTURE | Population assignment analyses                                  |
| PCA       | Principal component analyses                                    |
| FST       | Population differentiation analyses                             |

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
* BWA / BWA-MEM2
* SAMtools
* BCFtools
* VCFtools
* BEDTools
* ADMIXTURE
* GEMMA
* MSMC2
* Relate
* TreeMix
* Beagle
* IQ-TREE
* PLINK
* FastQC
* fastp
* Trimmomatic
* NVIDIA Parabricks
* R (various packages)

Please refer to each module directory for specific software requirements and installation instructions.

---

## Citation

If you use these pipelines in your research, please cite:

Chien, C.-C., Seiko, T., Muto, C., Ariga, H., Wang, Y.-C., Chang, C.-H., Sakai, H., Naito, K., & Lee, C.-R. (2025). *A single domestication origin of adzuki bean in Japan and the evolution of domestication genes*. **Science, 388**(6750), eads2871. https://doi.org/10.1126/science.ads2871

---

## Author

**Chih-Cheng Chien**

National Agriculture and Food Research Organization (NARO), Japan

National Taiwan University (NTU), Taiwan

### Research Interests

* Bioinformatics
* Computational genomics
* Population genomics
* Evolutionary genomics
* Genome editing
* Functional genomics
* HPC workflow automation
* Crop domestication and evolution
