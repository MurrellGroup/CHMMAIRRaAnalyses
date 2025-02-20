using CSV
using DataFrames
using Base.Filesystem

"""
File that inits IgDiscover runs for PCR_conditions libraries.
This assumes you've already run the genotyping step and have the databases in the PCR_conditions/IgM_genotyping_runs/ directory.
"""

# path to conda binary for running igdiscover
condabin = "/opt/miniconda3/condabin/conda"
# contains libraries dir with library folders 
# contains IgM_genotyping_runs with IgDiscover genotyping runs 
PCR_conditions_libraries_dir = "../data/libraries/PCR_conditions/libraries/"
PCR_conditions_genotyping_dir = "../data/igdiscover22/PCR_conditions_genotyping_runs/"
PCR_conditions_run_dir = "../data/igdiscover22/PCR_conditions/"

# aux file provided with repo
aux_path = "../data/human_gl.aux"
# these are the libraries made using our standard protocol
# we generate personalized genotypes using these and run assignments for all librarie using those genotypes
genotyping_libraries = ["D1-C1-IGM", "D2-C1-IGM", "D3-C1-IGM"]

donor2db_dir = Dict([(split(lib, "-")[1], joinpath(PCR_conditions_genotyping_dir, lib, "final", "database")) for lib in genotyping_libraries])

# relative path file finder
function recursive_file_finder(folder, substr)
	return vcat(filter(x->x != [], [[joinpath(el[1], file) for file in el[3] if occursin(substr, file)] for el in walkdir(folder)])...)
end

for donor in collect(keys(donor2db_dir))
	db_dir = donor2db_dir[donor]
	R1_files = recursive_file_finder(joinpath(PCR_conditions_libraries_dir, donor), "_R1_001.fastq.gz")
	for R1_file in R1_files
		run_dir = joinpath(PCR_conditions_run_dir, split(R1_file, "/")[end-1])
		init_cmd = `$(condabin) run -n igdiscover igdiscover init --database $(db_dir) --reads1 $(R1_file) $(run_dir)`
		println(init_cmd)
		run(init_cmd)
		cp("igdiscover_yamls/Chernyshev_2025/PCR_conditions_igdiscover.yaml", joinpath(run_dir, "igdiscover.yaml"), force = true)
		cp(joinpath(rsplit(db_dir, "/", limit = 3)[1], "database", "human_gl.aux"), joinpath(run_dir, "database", "human_gl.aux"))
	end
end
