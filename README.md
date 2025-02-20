# Analysis notebooks and scripts for "Detection of PCR chimeras in adaptive immune receptor repertoire sequencing using hidden Markov models"

## Description

Notebooks and scripts used to produces the figures and tables in the CHMMAIRRa paper. The figures generated from simulated data require only the databases already included in the repository. The figures generated from real data require preprocessing the datasets with [IgDiscover](https://gitlab.com/gkhlab/igdiscover22) according to the instructions in the [IgDiscover_preprocessing](./IgDiscover_preprocessing/) folder.

## Dependencies

Julia 1.10.5 for running the notebooks.
 - I recommend using [juliaup](https://github.com/JuliaLang/juliaup) to install Julia

 - All Julia package dependencies are listed in the [Manifest.toml](Manifest.toml) file.

 - Re-create the environment with the following commands in julia:
	```
	using Pkg; Pkg.activate("."); Pkg.instantiate()
	```

[IgDiscover v1.0.4](https://gitlab.com/gkhlab/igdiscover22) for preprocessing the real datasets.

## Notebooks
- [PCR_conditions.ipynb](./notebooks/PCR_conditions.ipynb) : Plots this paper's PCR parameter modification dataset.
- [databases.ipynb](./notebooks/databases.ipynb) : Plots pairwise edit distances between database V alleles.
- [simulations.ipynb](./notebooks/simulations.ipynb) : Simulation of TRB and IGH VDJ datasets to produce ROCs.
- [recombination.ipynb](./notebooks/recombination.ipynb) : Plots recombination information from real datasets. Produces the heatmaps, recombination percentage scatterplots, and database subsampling scatterplots.
- [benchmark_speed.ipynb](./notebooks/benchmark_speed.ipynb) : Benchmarks the speed of CHMMAIRRa, USEARCH, and VSEARCH on real and simulated datasets.
- [lineages.ipynb](./notebooks/lineages.ipynb) : Plots lineage information from a real dataset.
- [summarize_seqcounts.ipynb](./notebooks/summarize_seqcounts.ipynb) : Gathers sequence count data from all datasets.

## Scripts
- [run_CHMMAIRRa.sh](./scripts/run_CHMMAIRRa.sh) : Runs CHMMAIRRa on all 5 real datasets in the paper (4 published and 1 new).
- [run_CHMMAIRRa_db_subsampling.jl](./scripts/run_CHMMAIRRa_db_subsampling.jl) : Runs CHMMAIRRa on specific real TCR and IGH libraries with subsampled databases.
- [run_benckmarks.jl](./scripts/run_benckmarks.jl) : Runs CHMMAIRRa and uchime on varying sizes of simulated IGH and TRB datasets.

## IgDiscover preprocessing
- [IgDiscover_preprocessing](./IgDiscover_preprocessing/) : This folder contains instructions for preprocessing the real datasets with IgDiscover. One .md file for each of 5 datasets. Also contains descriptions of where to find the raw fastq data.