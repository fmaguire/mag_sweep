# MAG Sweep

Recover MAGs using a parameter sweep of settings implemented as a nextflow workflow.

Each tool can be entered into `run_params.csv` with each of the settings wanted and label on each line.
This is a combinatorial approach so for every assembler and parameter setting every binner and binner setting will be used to generate MAGs.

Currently implements the following metagenome assemblers:

- `metaSPAdes`
- `idba_ud`
- `megahit`

And the following MAG binners:

- `metabat2`
- `maxbin2`
- `concoct`

## Installation

Requires a working `conda` and `nextflow` (see respective projects for specific dependencies)

    curl -s https://get.nextflow.io | bash

All specific tool dependencies are installed via conda using the envs in `conda_envs/`

## Running

- Edit `nextflow.config` to select number of cores and input files to use
- Edit `run_params.csv` to add any parameters you want to add for the supported
tools above
- Execute the workflow (add `-resume` flag to rerun just incomplete parts of failed runs):

    ./nextflow run mag_sweep.nf

- Output can be found in `results/` by default

    results
    ├── 0.trimmed_reads
    ├── 1a_metaspades
    ├── 1b_idba_ud
    ├── 1c_megahit
    ├── 2_coverages
    ├── 3a_metabat2
    │   ├── idba_ud_default_metabat2_default
    │   ├── idba_ud_default_metabat2_low_score
    │   ├── idba_ud_default_metabat2_more_edges
    │   ├── idba_ud_default_metabat2_small_bins
    │   ├── megahit_default_metabat2_default
    │   ├── megahit_default_metabat2_low_score
    │   ├── megahit_default_metabat2_more_edges
    │   ├── megahit_default_metabat2_small_bins
    │   ├── megahit_large_metabat2_default
    │   ├── megahit_large_metabat2_low_score
    │   ├── megahit_large_metabat2_more_edges
    │   ├── megahit_large_metabat2_small_bins
    │   ├── megahit_sensitive_metabat2_default
    │   ├── megahit_sensitive_metabat2_low_score
    │   ├── megahit_sensitive_metabat2_more_edges
    │   ├── megahit_sensitive_metabat2_small_bins
    │   ├── metaspades_default_metabat2_default
    │   ├── metaspades_default_metabat2_low_score
    │   ├── metaspades_default_metabat2_more_edges
    │   ├── metaspades_default_metabat2_small_bins
    │   ├── metaspades_plasmidspades_metabat2_default
    │   ├── metaspades_plasmidspades_metabat2_low_score
    │   ├── metaspades_plasmidspades_metabat2_more_edges
    │   └── metaspades_plasmidspades_metabat2_small_bins
    ├── 3b_maxbin2
    │   ├── idba_ud_default_maxbin2_default
    │   ├── idba_ud_default_maxbin2_high_iteration
    │   ├── idba_ud_default_maxbin2_low_prob
    │   ├── megahit_default_maxbin2_default
    │   ├── megahit_default_maxbin2_high_iteration
    │   ├── megahit_default_maxbin2_low_prob
    │   ├── megahit_large_maxbin2_default
    │   ├── megahit_large_maxbin2_high_iteration
    │   ├── megahit_large_maxbin2_low_prob
    │   ├── megahit_sensitive_maxbin2_default
    │   ├── megahit_sensitive_maxbin2_high_iteration
    │   ├── megahit_sensitive_maxbin2_low_prob
    │   ├── metaspades_default_maxbin2_default
    │   ├── metaspades_default_maxbin2_high_iteration
    │   ├── metaspades_default_maxbin2_low_prob
    │   ├── metaspades_plasmidspades_maxbin2_default
    │   ├── metaspades_plasmidspades_maxbin2_high_iteration
    │   └── metaspades_plasmidspades_maxbin2_low_prob
    └── 3c_concoct
        ├── idba_ud_default_concoct_big_kmer
        ├── idba_ud_default_concoct_default
        ├── idba_ud_default_concoct_low_clusters
        ├── megahit_default_concoct_big_kmer
        ├── megahit_default_concoct_default
        ├── megahit_default_concoct_low_clusters
        ├── megahit_large_concoct_big_kmer
        ├── megahit_large_concoct_default
        ├── megahit_large_concoct_low_clusters
        ├── megahit_sensitive_concoct_big_kmer
        ├── megahit_sensitive_concoct_default
        ├── megahit_sensitive_concoct_low_clusters
        ├── metaspades_default_concoct_big_kmer
        ├── metaspades_default_concoct_default
        ├── metaspades_default_concoct_low_clusters
        ├── metaspades_plasmidspades_concoct_big_kmer
        ├── metaspades_plasmidspades_concoct_default
        └── metaspades_plasmidspades_concoct_low_clusters
