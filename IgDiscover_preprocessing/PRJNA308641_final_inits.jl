using CSV
using DataFrames
using Base.Filesystem

"""
File that inits IgDiscover assignment-only (0 iterations) runs for PRJNA308641 IgM and IgG libraries.
IgDiscover only provides accurate genotypes for IgM libraries, 
so we figure out which donor each IgG library belongs to and use the union of all its IgM library IgDiscover genotypes.
NOTE:
Based on jaccard indices of V genotypes in SRR3099483 and SRR3099450, the donor labels on the two libraries were switched!
SRR3099450 is actually from donor 1848 and SRR3099483 is from donor 1380. I correct this in the script below. 
All other libraries labeled with the same donor have high V genotype jaccard indices and therefore I think are labelled correctly. 
"""

sra_info_path = "PRJNA308641_SraRunInfo.csv"
condabin = "/home/mchernys/miniforge3/condabin/conda"
aux_path = "../data/human_gl.aux"
PRJNA308641_genotyping_dir = "../data/igdiscover22/PRJNA308641_genotyping_runs/IgM/"
PRJNA308641_assignment_dir = "../data/igdiscover22/PRJNA308641/"
PRJNA308641_libraries_dir = "../data/libraries/PRJNA308641/"

# Read and process SRA info
sra_info = CSV.read(sra_info_path, DataFrame)
sra_info.donor = map(x -> split(x, "_")[2], sra_info.LibraryName)
sra_info.chain = map(x -> split(x, "_")[end], sra_info.LibraryName)

# Correct donor labels
sra_info[sra_info.Run .== "SRR3099483", :donor] .= "1380"
sra_info[sra_info.Run .== "SRR3099450", :donor] .= "1848"

# Create donor databases
for donor_group in groupby(sra_info, :donor)
    donor = donor_group.donor[1]
    println(donor)
    db_dir = joinpath(PRJNA308641_genotyping_dir, "databases", donor)
    if !isdir(db_dir)
        mkpath(db_dir)
    end
    
    # Get IgM database paths
    igm_runs = donor_group[donor_group.chain .== "IgM", :Run]
    db_dirs = [joinpath(PRJNA308641_genotyping_dir, run, "final", "database") for run in igm_runs]
    V_dbs = [joinpath(db_dir, "V.fasta") for db_dir in db_dirs]
    D_dbs = [joinpath(db_dir, "D.fasta") for db_dir in db_dirs]
    J_dbs = [joinpath(db_dir, "J.fasta") for db_dir in db_dirs]
    
    # Run union commands
    for (gene, dbs) in [("V", V_dbs), ("D", D_dbs), ("J", J_dbs)]
        union_cmd = [condabin, "run", "-n", "igdiscover", "igdiscover", "union"]
        append!(union_cmd, dbs)
        union_cmd = `$union_cmd`
        output = pipeline(union_cmd, joinpath(db_dir, "$gene.fasta"))
        run(output)
    end
    cp(aux_path, joinpath(db_dir, "human_gl.aux"), force = true)
end

mkpath(PRJNA308641_assignment_dir)
mkpath(joinpath(PRJNA308641_assignment_dir, "IgM"))
mkpath(joinpath(PRJNA308641_assignment_dir, "IgG"))

# Initialize IgDiscover runs
for donor_group in groupby(sra_info, :donor)
    donor = donor_group.donor[1]
    db_dir = joinpath(PRJNA308641_genotyping_dir, "databases", donor)
    for row in eachrow(donor_group)
        init_cmd = `$condabin run -n igdiscover igdiscover init --database $db_dir --reads1 $(PRJNA308641_libraries_dir)/$(row.chain)/trimmed/$(row.Run)/cutadapt-$(row.Run)_1.fastq.gz $(PRJNA308641_assignment_dir)/$(row.chain)/$(row.Run)`
        run(init_cmd)
        cp("igdiscover_yamls/Galson_2015/PRJNA308641_igdiscover.yaml", joinpath(PRJNA308641_assignment_dir, "$(row.chain)", "$(row.Run)", "igdiscover.yaml"), force = true)
    end
end
