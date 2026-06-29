# GATK Pipeline v4.0

Automated and reproducible whole-genome sequencing (WGS) variant calling pipeline for High-Performance Computing (HPC) environments.

This pipeline provides end-to-end automation for NGS data processing, variant discovery, joint genotyping, variant filtering, and downstream analyses. It supports both CPU- and GPU-based workflows and is optimized for large-scale genomic studies.

---

## Features

* Automated end-to-end WGS analysis
* One-command execution
* HPC support (PBS, PBS Pro, and Slurm)
* Automatic job submission and dependency management
* GPU acceleration using NVIDIA Parabricks
* Incremental sample integration
* Germline CNV analysis
* WES and RAD-seq support
* Joint genotyping using CombineGVCFs or GenomicsDBImport
* Variant Quality Score Recalibration (VQSR)
* Scalable to thousands of samples

---

## Supported Workflows

* GATK4 pipeline (**recommended**)
* GPU-accelerated GATK4 pipeline (Parabricks)
* Legacy GATK3 pipeline
* Germline CNV analysis
* Whole-genome sequencing (WGS)
* Whole-exome sequencing (WES)
* RAD-seq / DArT-seq data

---

## Workflow Overview

```text
Step 0  Reference indexing
Step 1  FASTQ QC and trimming
Step 2  Read mapping and GVCF generation
Step 3  GVCF combination / GenomicsDBImport
Step 4  Joint genotyping
Step 5  Variant filtering
Step 6  Bi-allelic SNP selection and additional filtering
Step 7  Site depth calculation
```

Typical workflow:

```text
FASTQ
 ↓
QC / Trimming
 ↓
Alignment
 ↓
Duplicate removal
 ↓
Variant calling (GVCF)
 ↓
Joint genotyping
 ↓
Variant filtering
 ↓
Final VCF
```

---

## Requirements

### Required Software

* bwa or bwa-mem2
* samtools
* bcftools
* vcftools
* fastp
* GATK4

### Optional Software

* NVIDIA Parabricks
* Trimmomatic
* VEP
* R
* ggplot2

### Additional Requirements for CNV Mode

* Recent GATK4 installation with support for:

  * DetermineGermlineContigPloidy
  * GermlineCNVCaller
* VEP
* R (optional, for plotting)

---

## Environment Setup

Please ensure that all required software is installed and accessible through your `$PATH`.

Example:

```bash
export PATH=/path/to/software/bwa:$PATH
export PATH=/path/to/software/samtools:$PATH
export PATH=/path/to/software/bcftools:$PATH
export PATH=/path/to/software/vcftools/bin:$PATH
```

Recommended additions to `~/.bashrc`:

```bash
export PERL5LIB=/path/to/software/vcftools/perl/
export BCFTOOLS_PLUGINS=/path/to/software/bcftools/plugins/
export TRIMMO=/path/to/trimmomatic.jar
export LC_CTYPE="en_US.UTF-8"

export gatk4=/path/to/GATK_X.X.X/gatk
alias gatk4='/path/to/GATK_X.X.X/gatk'
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/aliveben2k/GATK_pipeline.git
```

Please make sure that:

* `qsub_subroutine.pl` is located in `$HOME`
* `filter_vcf_4.0.pl` is located in the same directory as `GATK_pipeline.pl`

---

## Quick Start

### Standard GATK4 Workflow

```bash
perl GATK_pipeline.pl \
    -p FASTQ_FOLDER \
    -r reference.fa \
    -g GROUP_NAME \
    -exc
```

### GPU Workflow (Parabricks)

```bash
perl GATK_pipeline.pl \
    -p FASTQ_FOLDER \
    -r reference.fa \
    -g GROUP_NAME \
    -gpu \
    -exc
```

### Legacy GATK3 Workflow

```bash
perl GATK_pipeline.pl gatk3 \
    -p FASTQ_FOLDER \
    -r reference.fa \
    -g GROUP_NAME \
    -exc
```

---

## Advanced Usage

### Run a Single Step

```bash
perl GATK_pipeline.pl \
    -p PATH \
    -r reference.fa \
    -sp 3s \
    -exc
```

### Run Multiple Steps

```bash
perl GATK_pipeline.pl \
    -p PATH \
    -r reference.fa \
    -g GROUP_NAME \
    -sp 1p \
    -esp 4 \
    -exc
```

### Dry Run Mode

Generate job scripts without submitting jobs:

```bash
perl GATK_pipeline.pl \
    -p PATH \
    -r reference.fa
```

### Background Execution (Recommended)

```bash
nohup perl GATK_pipeline.pl \
    -p PATH \
    -r reference.fa \
    -g GROUP_NAME \
    -exc > Log.txt 2>&1 &
```

Monitor progress:

```bash
tail -f Log.txt
```

---

## CNV Workflow

### Cohort Mode

```bash
perl GATK_pipeline_v4.0_gpu.pl \
    -p PATH \
    -r reference.fa \
    -g GROUP_NAME \
    -cnv \
    -pps PLOIDY_PRIORS.tsv \
    -sp 2p \
    -esp 3 \
    -exc
```

### Case Mode

```bash
perl GATK_pipeline_v4.0_gpu.pl \
    -p PATH \
    -r reference.fa \
    -g GROUP_NAME \
    -cnv \
    -pps PLOIDY_PRIORS.tsv \
    -chf COHORT_MODEL_FOLDER \
    -sp 2p \
    -esp 3 \
    -exc
```

---

## Pipeline Steps

| Step | Description                          |
| ---- | ------------------------------------ |
| 0    | Reference indexing                   |
| 1    | FASTQ quality control and trimming   |
| 2    | Read mapping and GVCF generation     |
| 3    | GVCF combination or GenomicsDBImport |
| 4    | Joint genotyping                     |
| 5    | Variant filtering                    |
| 6    | SNP and genotype filtering           |
| 7    | Site depth calculation               |

---

## Important Notes

* **GATK4 is strongly recommended.**
* GATK3-generated GVCFs should not be used in GATK4 workflows.
* Use `-d` (GenomicsDBImport) when sample size exceeds ~1500.
* If the reference contains many contigs, use `-ns`.
* The pipeline automatically checks existing files and can resume interrupted analyses.
* Use `Check_log.pl` to identify incomplete jobs or failed outputs.

---

## Documentation

Detailed descriptions of command-line arguments are available in:

```text
docs/arguments.md
```

---

## Citation

If you use this pipeline in your research, please cite:

Chien, C.-C. et al.

---

## Author

**Chih-Cheng Chien**

National Agriculture and Food Research Organization (NARO), Japan

Research interests:

* Bioinformatics
* Computational genomics
* Population genomics
* Evolutionary genomics
* Genome editing
* HPC workflow automation
