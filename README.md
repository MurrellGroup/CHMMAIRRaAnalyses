# Analysis notebooks and scripts for "Detection of PCR chimeras in adaptive immune receptor repertoire sequences"

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

[USEARCH v11.0.667_i86linux32](drive5.com/usearch/) comparison method.

[VSEARCH v2.29.1_linux_x86_64](github.com/torognes/vsearch) comparison method.

[IgDiscover v1.0.4](https://gitlab.com/gkhlab/igdiscover22) for preprocessing the real datasets.

[MAFFT v7.490](mafft.cbrc.jp/alignment/software/) for reference database alignment.

[Muscle v5.3](https://anaconda.org/bioconda/muscle) for reference database alignment.

[ART 2016-06-05](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) to simulate MiSeq noise in TRB sequences.

[Shazam v1.2.0](https://www.rdocumentation.org/packages/shazam/versions/1.2.0) to simulate SHM in IGH sequences.

## Simulated data evaluation (Fig. 3, Fig. S1)

These scripts requires the Julia environment, ART, Shazam, VSEARCH, and USEARCH. 

- [simulate_IGH_shazam.jl](./scripts/simulate_IGH_shazam.jl) : Simulates IGH V, D, and J datasets with SHM added by shazam's shmulateSeq.
- [simulate_TRB_art_illumina.jl](./scripts/simulate_TRB_art_illumina.jl) : Simulates TRB V, D, and J datasets with sequencing noise added by art_illumina from ART.
- [ROCs.ipynb](./notebooks/ROCs.ipynb) : Generates ROCs for simulated TRB and IGH V, D, and J datasets.
- [run_benckmarks.jl](./scripts/run_benckmarks.jl) : Runs CHMMAIRRa and uchime on varying sizes of simulated IGH and TRB datasets.
- [benchmark_speed.ipynb](./notebooks/benchmark_speed.ipynb) : Plots the speed of CHMMAIRRa, USEARCH, and VSEARCH on simulated (and real) datasets.

## Main real data analysis

The analysis of the real data involves running IgDiscover on 319 libraries with dataset-specific settings, so this takes some doing. Requires IgDiscover and the Julia environment.
- [IgDiscover_preprocessing](./IgDiscover_preprocessing/) : This folder contains instructions for preprocessing the real datasets with IgDiscover. One .md file for each of 5 datasets. Also contains descriptions of where to find the raw fastq data.
- [run_CHMMAIRRa.jl](./scripts/run_CHMMAIRRa.jl) : Runs CHMMAIRRa on all 5 real datasets in the paper (4 published and 1 new).
- [run_CHMMAIRRa_db_subsampling.jl](./scripts/run_CHMMAIRRa_db_subsampling.jl) : Runs CHMMAIRRa on specific real TCR and IGH libraries with subsampled databases (for Fig. 4).
- [recombinations.ipynb](./notebooks/recombinations.ipynb) : Plots recombination information from real datasets. Produces the heatmaps (Fig. 7), recombination percentage scatterplots (Fig 6), and database subsampling scatterplots (Fig. 4) as well as Supplementary Figures S2, S3, and S7.
- [PCR_conditions.ipynb](./notebooks/PCR_conditions.ipynb) : Plots this paper's PCR parameter modification dataset (Fig. 5).
- [lineages.ipynb](./notebooks/lineages.ipynb) : Plots lineage information from a real dataset (Fig. 1).
- [summarize_seqcounts.ipynb](./notebooks/summarize_seqcounts.ipynb) : Gathers sequence count data from all datasets (Supplementary Data 1).

## Other analyses
- [run_CHMMAIRRa_varying_V_DFR.jl](./scripts/run_CHMMAIRRa_varying_V_DFR.jl) : Runs CHMMAIRRa at varied minimum differences from reference (DFR) settings (for Supplementary Figure S2). 
- [compare_alignment_settings.jl](./scripts/compare_alignment_settings.jl) : Runs a set of libraries with varying V database alignment methods (for Supplementary Figure S3).
- [databases.ipynb](./notebooks/databases.ipynb) : Plots pairwise edit distances between database V alleles (Supplementary Figure S6).
- [run_CHMMAIRRa_Js.jl](./scripts/run_CHMMAIRRa_Js.jl) : Runs CHMMAIRRa on J segments for a few libraries to get a grasp on their frequency.