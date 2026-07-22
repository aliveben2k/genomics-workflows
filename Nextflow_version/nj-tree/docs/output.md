# Output

For `--output_prefix adzuki`, the final output directory contains:

```text
adzuki.mtx                 Final labeled distance matrix
adzuki_final.npz           NumPy archive of merged matrices
adzuki_pcoa.txt            PCoA coordinates
adzuki_pcoa.pdf            PCoA scatter plot
adzuki_tree.newick         Neighbor-joining tree in Newick format
adzuki_tree.nex            Neighbor-joining tree in Nexus format
adzuki_tree_pcoa.pkl       Python PCoA/tree result bundle
```

Intermediate outputs are retained for reproducibility and restart inspection:

```text
intermediate/tables/       Compressed genotype-table chunks and manifest
intermediate/matrices/     Per-chunk available-site/difference matrices
pipeline_info/             Nextflow trace, timeline, report, and DAG
```
