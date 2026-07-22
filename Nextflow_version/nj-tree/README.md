# aliveben2k/genomics-workflows (Nextflow_version/nj-tree)

An nf-core-template-based Nextflow DSL2 pipeline for chunked pairwise genomic
distance calculation, principal coordinates analysis (PCoA), and
neighbor-joining tree reconstruction.

## Workflow

1. Convert a diploid, biallelic VCF/BCF into compressed genotype-table chunks.
2. Scatter one independent pairwise-distance process per chunk.
3. Gather the available-site and difference matrices.
4. Calculate the final normalized distance matrix.
5. Generate PCoA coordinates/plot and Newick/Nexus NJ trees.

Nextflow automatically submits the scattered pairwise processes to the selected
executor and waits for all successful outputs before running finalization.

## Quick start

```bash
nextflow run . \
    --input variants.vcf.gz \
    --pop populations.tsv \
    --outdir results \
    --output_prefix adzuki \
    --window_size 10000 \
    -profile conda,local
```

Slurm:

```bash
nextflow run . \
    --input variants.vcf.gz \
    --pop populations.tsv \
    --outdir results \
    --output_prefix adzuki \
    -profile conda,slurm \
    -resume
```

PBS Pro:

```bash
nextflow run . \
    --input variants.vcf.gz \
    --outdir results \
    --output_prefix adzuki \
    -profile conda,pbspro \
    -resume
```

The population file is optional. Without `--pop`, output distances are
calculated among individual samples.

See [usage documentation](docs/usage.md) and
[output documentation](docs/output.md).

## Requirements

- Nextflow 25.10.4 or newer
- Conda/Mamba when using the bundled `conda` or `mamba` profiles, or an
  environment that already provides the packages in `environment.yml`

This project was created with the official nf-core pipeline template and uses
local nf-core-style modules for the project-specific Python programs.
