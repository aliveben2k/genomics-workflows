# Usage

## Required parameters

```text
--input           Input .vcf, .vcf.gz, or .bcf
--outdir          Output directory
```

## Optional parameters

```text
--pop             Two-column tab-separated sample/population file
--output_prefix   Output filename prefix (default: nj-tree)
--window_size     Retained variants per chunk (default: 10000)
```

The population file may contain a header beginning with `ind`, `id`, `taxa`,
or `samp`. Column one is the VCF sample ID and column two is the population.

## Local execution

```bash
nextflow run . \
    --input variants.vcf.gz \
    --outdir results \
    -profile conda,local
```

## Slurm

```bash
nextflow run . \
    --input variants.vcf.gz \
    --pop populations.tsv \
    --outdir results \
    -profile conda,slurm \
    -resume
```

## PBS and PBS Pro

```bash
# PBS/Torque
nextflow run . --input variants.vcf.gz --outdir results -profile conda,pbs

# PBS Pro
nextflow run . --input variants.vcf.gz --outdir results -profile conda,pbspro
```

Cluster-specific queues, projects, resource limits, and module-loading commands
should be provided in an institutional Nextflow config:

```bash
nextflow run . -c cluster.config -profile conda,slurm ...
```

Pairwise chunks are independent Nextflow tasks. Their maximum concurrency is
controlled by the executor and may be limited in a custom config with
`process.queueSize`.

Use `-resume` when restarting so successfully completed chunks are reused.
