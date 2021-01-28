# MAG Sweep

Recover MAGs using a parameter sweep of settings implemented as a nextflow workflow.

Each tool can be entered into `run_params.csv` with each of the settings wanted and label on each line.
This is a combinatorial approach so for every assembler and parameter setting every binner and binner setting will be used to generate MAGs.

