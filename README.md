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


