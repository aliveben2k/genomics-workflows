# nf-core/gatkpopulation: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Read trimming](#read-trimming) - fastp or Trimmomatic output
- [FastQC](#fastqc) - Trimmed-read QC
- [Reference preparation](#reference-preparation) - FASTA, dictionary, and BWA indexes
- [Read preprocessing](#read-preprocessing) - CPU BWA/GATK or GPU Parabricks output
- [Per-sample GVCFs](#per-sample-gvcfs) - Indexed GVCFs from GATK4 or Parabricks
- [GenomicsDB](#genomicsdb) - One all-sample workspace per reference contig
- [Subset GVCFs](#subset-gvcfs) - Optional selected-sample GVCF per contig
- [Raw variant calls](#raw-variant-calls) - Per-contig multi-sample VCFs
- [Filtered variant calls](#filtered-variant-calls) - Configurable bcftools view output
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Read trimming

- `01-fq_trim/`: trimmed reads and fastp or Trimmomatic reports.

Retained reads use `${sample}.fastq.trimmed.gz` for single-end data or
`${sample}_R1.fastq.trimmed.gz` and `${sample}_R2.fastq.trimmed.gz` for
paired-end data.

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### Reference preparation

- `00-ref_gatk4/`: FASTA `.fai`, GATK sequence dictionary, and the traditional
  BWA index generated when `--bwa_index` is not supplied.

### Read preprocessing

- `02-bam/`: CPU BWA-MEM or Parabricks alignment, BAM index, optional duplicate
  metrics and BQSR table, final BAM, and the `samtools stats` report.

When no known-sites VCF is supplied, the final output is the indexed BAM and
BQSR is skipped.

WGS and WES expose `${sample}_aln_sort_MD.bam` after duplicate marking. RAD,
GBS, ddRAD, and amplicon modes skip duplicate marking and expose
`${sample}_aln_sort.bam`. ApplyBQSR uses a `.tmp.bam` internally before
restoring the applicable final name. Statistics follow the same basename:
`${sample}_aln_sort_MD.stats` or `${sample}_aln_sort.stats`.

### Per-sample GVCFs

- `03-sample_gvcf/`: `${sample}.g.vcf.gz` and
  `${sample}.g.vcf.gz.tbi` for each final BAM.

CPU runs use GATK4 HaplotypeCaller with GVCF reference-confidence output. GPU
runs use Parabricks `haplotypecaller --gvcf` followed by `pbrun indexgvcf`.
When GPU fallback is enabled, an unsuccessful GPU sample is called with GATK4.

### GenomicsDB

- `04-genomicsDB/`: one `${contig}.genomicsdb/` workspace per contig in
  the reference FASTA index.

Each workspace contains all indexed sample GVCFs for its contig and is created
by an independent GATK4 `GenomicsDBImport` job. The generated sample-name map is
an internal workflow file.

### Subset GVCFs

- `05.1-subset/`: `${contig}.subset.g.vcf.gz` and
  `${contig}.subset.g.vcf.gz.tbi`.

These files are created only when `--subset_samples` is supplied. The input CSV
must contain exactly one column named `sample`. GATK4 `SelectVariants` reads
each GenomicsDB workspace and creates a flat GVCF containing only those
samples. This operation does not modify the source GenomicsDB.

### Raw variant calls

- `05-raw_vcf/`: per-contig multi-sample variant calls.
- `05-raw_vcf/all.raw.vcf.gz`: optional gathered raw VCF and `.tbi`.

All profiles use GATK4 `GenotypeGVCFs` and write `${contig}.raw.vcf.gz` plus
`${contig}.raw.vcf.gz.tbi`. By default, GenotypeGVCFs reads each per-contig
GenomicsDB directly. With `--subset_samples`, it instead reads the corresponding
`05.1-subset` GVCF. The `gpu` profile is used only for earlier preprocessing
and per-sample GVCF generation.

With `--all_sites`, GATK4 includes non-variant loci and writes
`${contig}.all_sites.raw.vcf.gz` plus `.tbi`. All-sites VCFs can be
substantially larger than regular variant-only VCFs.

With `--gather_vcfs 5`, GATK4 `GatherVcfs` combines the per-contig raw files
into `all.raw.vcf.gz` in reference contig order. Stage 6 then filters only this
gathered VCF. The published per-contig raw files remain available.

### Filtered variant calls

- `06-filtered_vcf/`: `${contig}.filtered.vcf.gz` and
  `${contig}.filtered.vcf.gz.tbi`.
- `06-filtered_vcf/all.filtered.vcf.gz`: filtered genome-wide VCF and `.tbi`
  when gathering is enabled.

The default `--filter_vcf biallele` preset retains biallelic SNPs with less
than 10% missing data, `QUAL >= 30`, and `INFO/DP > 3`:

```text
-v snps -m2 -M2 -i 'F_MISSING < 0.1 && QUAL >= 30 && INFO/DP > 3'
```

The `--filter_vcf monobi` preset retains monomorphic reference sites and
biallelic SNPs, removes indels, and applies `QUAL >= 30` only to SNP records:

```text
-m1 -M2 -i 'F_MISSING < 0.1 && INFO/DP > 3 && (TYPE="ref" || (TYPE="snp" && QUAL >= 30))'
```

Use `monobi` with `--all_sites` unless the stage-6 input already contains
monomorphic records.

With `--gather_vcfs 6`, each contig is filtered independently before GATK4
`GatherVcfs` combines the results into `all.filtered.vcf.gz`. With
`--gather_vcfs 5`, `all.filtered.vcf.gz` is produced by filtering the gathered
`all.raw.vcf.gz`.

Use `--bcftools_view_args` or a Nextflow parameter file to override either
preset and access other `bcftools view` filtering or subsetting options. Output
remains bgzip-compressed and tabix-indexed.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
