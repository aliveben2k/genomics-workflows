# Relate_v8 Pipeline

Wrapper pipeline for running Relate from VCF input, with optional downstream support for:

- population size estimation
- DetectSelection / DPS
- CLUES2
- TreeViewSamples

This pipeline is tested with **Relate v1.2.4**.

## Official Relate Documentation

- [Relate documentation](https://myersgroup.github.io/relate/)

## Important Compatibility Note

Before using TreeViewSamples with this wrapper, replace the official Relate v1.2.4 files:

- `scripts/TreeView/TreeViewSample.sh`
- `scripts/TreeView/treeview_sample.R`

with the **modified versions provided in this pipeline/workspace**.

## Main Features

- accepts either:
  - a folder of `.vcf` / `.vcf.gz` files
  - a text file listing VCF paths
- compresses plain `.vcf` files automatically
- checks and creates `.tbi` indexes for `.vcf.gz` files
- splits a multi-contig VCF into per-contig input VCFs when needed
- uses `vcf2anc_vcf_threads.pl`
- uses the main `-t` option as the shared thread control
- uses **SHAPEIT5** instead of pseudo-phasing
- validates and reorders `*.poplabels` to match the Relate `.sample` file
- supports haploid mode with `-hap`
- supports removal of ancestral samples with `-nka`
- supports balanced repeated sampling with `-rs`
- supports repeated CLUES2 runs with `-rc`
- automatically summarizes CLUES2 and TreeViewSamples outputs
- uses CLUES2 summary outputs for AIC-based evaluation and broken-stick breakpoint analysis

## Required Inputs

### Standard VCF-based workflow

- VCF file(s): folder or list file
- poplabels file: supplied with `-pop`
- ancestral sample ID list: supplied with `-al`
- genetic map file: supplied with `-map`

### Optional inputs

- mask fasta: `-mask`
- existing anc/mut inputs: `-am`

If `-am` is used, the pipeline starts from existing `.anc/.mut` files, so `-vcf`, `-map`, and `-al` are not required.

## Poplabels Requirements

Header must be exactly:

```text
sample population group sex
```

Rules used by the pipeline:

- if `-hap` is used, `sex` is rewritten to `1`
- otherwise `sex` is rewritten to `NA`
- sample order is checked against the Relate `.sample` input
- missing samples in poplabels will stop the run with an error
- if `-nka` is used, ancestral samples are removed from poplabels before downstream steps
- if `-rm` is used, removed samples are also filtered from poplabels

## Phasing Behavior

The pipeline no longer uses pseudo-phasing in `vcf2anc_vcf_threads.pl`.

For VCF input:

1. an intermediate VCF is prepared
2. `bcftools +fill-tags` adds AC/AN tags
3. SHAPEIT5 is run through conda
4. SHAPEIT5 BCF output is converted into Relate `.haps` / `.sample`

The conda environment name currently used by the wrapper is:

```perl
$shapeit = 'shapeit5'
```

Example SHAPEIT5 command used by the workflow:

```bash
conda run -n shapeit5 SHAPEIT5_phase_common --input INPUT.vcf.gz --map CHR.gmap.map --region CHR --output PHASED.bcf --thread THREADS
```

If `-hap` is used, SHAPEIT5 still runs first, then the first phased allele is used to build one haplotype per sample.

## Basic Example

```bash
nohup perl Relate_v8.pl \
  -vcf VCF_FOLDER_OR_LIST \
  -pop POPLABELS_FILE \
  -map RECOMB_MAP_FILE \
  -al ANCESTOR_ID_LIST \
  -o OUTPUT_PATH \
  -t 8 \
  -exc > LOGFILE 2>&1 &
```

## Example Commands

### Standard diploid run

```bash
nohup perl Relate_v8.pl \
  -vcf vcfs.list \
  -pop samples.poplabels \
  -map rmap_all.txt \
  -al ancestral_ids.txt \
  -o Relate_out \
  -t 8 \
  -exc > Relate.log 2>&1 &
```

### Haploid-style run

```bash
nohup perl Relate_v8.pl \
  -vcf vcfs.list \
  -pop samples.poplabels \
  -map rmap_all.txt \
  -al ancestral_ids.txt \
  -hap \
  -o Relate_hap \
  -t 8 \
  -exc > Relate_hap.log 2>&1 &
```

### Run without keeping ancestral samples

```bash
nohup perl Relate_v8.pl \
  -vcf vcfs.list \
  -pop samples.poplabels \
  -map rmap_all.txt \
  -al ancestral_ids.txt \
  -nka \
  -o Relate_no_anc \
  -t 8 \
  -exc > Relate_no_anc.log 2>&1 &
```

### CLUES2 + TreeViewSamples together

```bash
nohup perl Relate_v8.pl \
  -vcf vcfs.list \
  -pop samples.poplabels \
  -map rmap_all.txt \
  -al ancestral_ids.txt \
  -clues \
  -tvs \
  -bp Chr01:100000-120000 \
  -o Relate_clues_tvs \
  -t 8 \
  -exc > Relate_clues_tvs.log 2>&1 &
```

### Balanced repeated sample sets

```bash
perl Relate_v8.pl \
  -vcf vcfs.list \
  -pop samples.poplabels \
  -map rmap_all.txt \
  -al ancestral_ids.txt \
  -clues \
  -bp Chr01:100000-120000 \
  -rs 10,100 \
  -rseed 123 \
  -t 8 \
  -exc
```

## Key Options

### VCF / sample handling

- `-t INT`: main thread count
- `-hap`: haploid output mode
- `-nka`: do not keep ancestral samples
- `-rm FILE`: remove listed samples
- `-mask FILE`: masked fasta input

### CLUES2 / TreeViewSamples

- `-clues`: prepare CLUES2 outputs
- `-tvs`: prepare TreeViewSamples outputs
- `-bp CHR:START-END[,CHR:START-END...]`: target region(s)
- `-rc INT`: repeat CLUES2 inference runs
- `-rs SAMPLE_NUM,REPEATS`: generate balanced repeated sample sets
- `-rseed INT`: random seed for `-rs`
- `-nat`: skip CLUES2 trajectory reconstruction / summary plotting
- `-tvs_debug`: save TreeViewSamples debug RDS bundles for figure redesign

### CLUES2 summary / rerun control

- `-rp clues`
- `-rp tvs`
- `-rp all`
- `-ow`
- `-npsi INT`: number of breakpoints for broken-stick analysis, default `1`

### Population size

- `-epsrp INT`: repeat EstimatePopulationSize MCMC runs
- default `-epsrp` is `0`
- `-bins LOWER,UPPER,STEPSIZE`: custom EPS bins

### CLUES2 plotting / inference

- `-cb FILE` or `-cb values`: custom CLUES time bins
- `-tco INT`: CLUES `tCutoff`
- `-d FLOAT`: CLUES dominance parameter

## Downstream Summary Behavior

### CLUES2

Unless `-nat` is used, the pipeline automatically generates summarized CLUES2 outputs, including:

- `AIC_all_runs.txt`
- `AIC_lowest.txt`
- `AIC_summary.txt`
- `final_all_loci.median.pdf`
- `final_all_loci.broken.stick.rda`
- `broken_stick_results.rda`
- `broken_stick_results.txt`

CLUES2 summary processing in this pipeline includes:

- AIC-based evaluation across repeat runs
- summarized allele-frequency trajectory plotting
- broken-stick model fitting on summarized CLUES2 trajectories
- breakpoint evaluation with `-npsi` controlling the number of breakpoints

The old `-rpp` behavior has been removed.

### TreeViewSamples

The pipeline automatically summarizes repeat outputs and writes merged TreeViewSamples summary files.

If `-tvs_debug` is used, extra debug RDS files are written containing:

- `summary_data`
- `output_prefix`

## Required Pipeline Files

- `Relate_v8.pl`
- `vcf2anc_vcf_threads.pl`
- `validate_poplabels_Relate.pl`
- `reorder_poplabels_by_sample_Relate.pl`
- `random_sets_Relate.pl` when `-rs` is used
- `convert_masked_fasta.pl` when `-mask` is used
- `recomb_spline_Relate.R`
- `plot_population_size_new.R`
- `DPS_plot.R` when DPS plotting is used
- `extract_derived_file.pl`
- `qsub_subroutine.pl` in `~/softwares/` or `~/`

## Local Helper Files Used by This Wrapper

### CLUES2 helper files

- `clues/RelateToCLUES.py`
- `clues/inference.py`
- `clues/all.freq.code.revision.clues2.R`
- `clues/cal_AICs_v2.R`
- `clues/broken_stick_analysis_v3_no_0.R`

### TreeView helper files

- `relate_v1.2.4/scripts/TreeView/TreeViewSample.sh`
- `relate_v1.2.4/scripts/TreeView/treeview_sample.R`

## External Software Requirements

- Relate v1.2.4
- bcftools with `+fill-tags`
- bgzip
- tabix
- conda environment `shapeit5`
- Python
- R / Rscript
- required R packages, including:
  - `segmented`
  - `ggplot2`
  - `cowplot`
  - `dplyr`
  - `reshape2`
  - `RColorBrewer`
- cluster submission environment when `-exc` is used

## Notes

- The main supported v8 workflow is `-pop` with optional `-hap`
- the pipeline renames contigs internally and writes `rename_chr.list`
- incomplete poplabels stop the run before Relate analysis
- use `--force` carefully for very large `-rc` values
