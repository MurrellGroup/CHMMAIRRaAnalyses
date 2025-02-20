"""
Script to run benchmarks for CHMMAIRRa and uchime on simulated and real datasets.

We simulate IGH and TRB genotypes from the OGRDB and KI TCR DB, generating 1000, 10000, and 100000 sequences for each.
We also benchmark CHMMAIRRa and uchime on TRB library PRJNA300878/TRB/1A and an IGG library SRR3099387.
"""


using Pkg 
Pkg.develop(path = "/home/mchernys/ben/chimera_detection/CHMMAIRRa.jl")
using CHMMAIRRa
Pkg.activate("..")
# using CHMMAIRRa
using Random, CSV, DataFrames, JSON
include("../src/simulate.jl")
include("../src/benchmark.jl")
include("../src/utils.jl")



# INPUTS
igdiscover_dir = "../../data/igdiscover22/"
CHMMAIRRa_bin = "../data/CHMMAIRRa_ubuntu-linux_x86-64_v0.0.1/bin/CHMMAIRRa" # obtained from https://github.com/MurrellGroup/CHMMAIRRa.jl/releases
# OUTPUTS
benchmark_dir = "../../outputs/timing_benchmarks/"



KITDB_dir = joinpath(@__DIR__, "../data/KI_TCR_DB_v0.0.1") # Provided with repository. See https://gkhlab.gitlab.io/tcr/sequences/
OGRDB_dir = joinpath(@__DIR__, "../data/OGRDB_human_IGH_9") # Provided with repository. See https://ogrdb.airr-community.org/germline_sets/Homo%20sapiens

datasets = ["TRB sim 1000", "TRB sim 10000", "TRB sim 100000", "TRB (PRJNA300878/TRB/1A)", "IGH sim 1000", "IGH sim 10000", "IGH sim 100000", "IGG (SRR3099387)"]
dataset2ind = Dict(zip(datasets, collect(1:length(datasets))));
random_seed = 888
rng = MersenneTwister(random_seed)

#### Create synthetic IGH genotype from OGRDB 

OGRDB_json = JSON.parsefile(joinpath(OGRDB_dir, "Homo_sapiens_IGH_VDJ_rev_9_ex.json"))
OGRDB_refnames, OGRDB_refseqs = read_fasta(joinpath(OGRDB_dir, "V.fasta"));

OGRDB_functional_alleles = [el["label"] for el in OGRDB_json["GermlineSet"][1]["allele_descriptions"] if el["functional"] & occursin("IGHV", el["label"])]
functional_inds = OGRDB_refnames .∈ Ref(OGRDB_functional_alleles)
OGRDB_genotype_refnames, OGRDB_genotype_refseqs = simulate_genotype(OGRDB_refnames[functional_inds], OGRDB_refseqs[functional_inds], rng)

#### Create synthetic TRB genotype from KI TCR DB

KITDB_TRB_refnames, KITDB_TRB_refseqs = read_fasta(joinpath(KITDB_dir, "TRB", "V.fasta"));
KITDB_TRB_genotype_refnames, KITDB_TRB_genotype_refseqs = simulate_genotype(KITDB_TRB_refnames, KITDB_TRB_refseqs, rng)


#### Benchmark CHMMAIRRa and uchime on varying sizes of simulated IGH and TRB datasets 

# set shm rates and query counts
simulated_genotypes = Dict("IGH" => (OGRDB_genotype_refnames, OGRDB_genotype_refseqs), "TRB" => (KITDB_TRB_genotype_refnames, KITDB_TRB_genotype_refseqs))
db_dirs = Dict("IGH" => OGRDB_dir, "TRB" => joinpath(KITDB_dir, "TRB"))
chain2receptor = Dict("IGH" => "IG", "TRB" => "TCR")
shm1,shm2 = 0.0,0.0
Ns = [1000, 10000, 100000]
sim_times = DataFrame()
chmmairra_is_compiled = false
for (chain, n) in Base.product(["TRB", "IGH"], Ns)
        mktempdir() do dir
        println("---------------------")
        println((chain, n))
        println("---------------------")

        # IgBLAST fails to align some chimeras, causing duplicate v_sequence_alignment values despite differing sequences. This is unfair for timing.
        # to solve this, I simulate a bit more sequences than needed, then subset the IgBLAST results based on unique v_sequence_alignments
        chimeric_names, chimeric_seqs, ref_names, ref_seqs, breakpoint_positions = random_chimeras(simulated_genotypes[chain][1], simulated_genotypes[chain][2], rng,
                                                                                                        min_pos = 1/4,
                                                                                                        max_pos = 3/4,
                                                                                                        n = n,
                                                                                                        min_shm1 = shm1,
                                                                                                        max_shm1 = shm1,
                                                                                                        min_shm2 = shm2,
                                                                                                        max_shm2 = shm2);
        # set up input paths
        db_fasta_path = joinpath(dir, "V.fasta")
        query_fasta_path = joinpath(dir, "queries_$(n).fasta")
        assignments_path = joinpath(dir, "assignments_$(n).tsv")
        augmented_assignments_path = joinpath(dir, "augmented_$(n).tsv")

        # write db 
        write_fasta(db_fasta_path, degap.(simulated_genotypes[chain][2]), seq_names = simulated_genotypes[chain][1])
        cp(joinpath(db_dirs[chain], "D.fasta"), joinpath(dir, "D.fasta"), force = true)
        cp(joinpath(db_dirs[chain], "J.fasta"), joinpath(dir, "J.fasta"), force = true)

        # write query for uchime
        write_fasta(query_fasta_path, degap.(chimeric_seqs), seq_names = chimeric_names)

        # write query dataframe for CHMMAIRRa
        run(pipeline(`conda run -n igdiscover igdiscover igblastwrap --threads $(Threads.nthreads()) $(dir) $(query_fasta_path)`, stdout = assignments_path))
        run(pipeline(`conda run -n igdiscover igdiscover augment $(dir) $(assignments_path)`, stdout = augmented_assignments_path))

        # make sure CHMMAIRRa is compiled before running benchmarks
        global chmmairra_is_compiled
        if ! chmmairra_is_compiled
            precompile_chmmairra(db_fasta_path, augmented_assignments_path)
            chmmairra_is_compiled = true
        end

        global chain2receptor    
        timing_tuple = time_all_methods(db_fasta_path, augmented_assignments_path, query_fasta_path,
                                        receptor = chain2receptor[chain],
                                        uchime_mode = "balanced",
                                        CHMMAIRRa_bin = CHMMAIRRa_bin)

        global sim_times
        sim_times = vcat(sim_times, DataFrame(
            "Method" => collect(keys(timing_tuple)),
            "Time (s)" => collect(values(timing_tuple)),
            "Dataset" => ["$(chain) sim $(n)" for i in 1:length(timing_tuple)]))

        println("Current timings $(timing_tuple)")
    end
end
CSV.write(joinpath(benchmark_dir, "sim_times.tsv"), sim_times, delim = "\t")

#### Benchmark CHMMAIRRa and uchime on real TRB and IGH data

igdiscover_dirs = [joinpath(igdiscover_dir, "PRJNA300878/TRB/1A/final"), joinpath(igdiscover_dir, "PRJNA308641/IgG/SRR3099387/final/")]
dataset_names = ["TRB (PRJNA300878/TRB/1A)", "IGG (SRR3099387)"]
chain2receptor = Dict("IGH" => "IG", "TRB" => "TCR")
chain = ["TRB", "IGH"]
real_data_times = DataFrame()
for (igdiscover_dir, chain, dataset_name) in zip(igdiscover_dirs, chain, dataset_names)
    println((igdiscover_dir, chain))
    query_tsv_path = joinpath(igdiscover_dir, "filtered.tsv")
    db_fasta_path = joinpath(igdiscover_dir, "database", "V.fasta")
    query_fasta_path = joinpath(igdiscover_dir, "queries.fasta")

    query_df = CSV.read(query_tsv_path, DataFrame, delim = "\t", select = ["sequence_id", "v_sequence_alignment"])
    write_fasta(query_fasta_path, degap.(convert(Vector{String}, query_df.v_sequence_alignment)), seq_names = query_df.sequence_id)

    timing_tuple = time_all_methods(db_fasta_path, query_tsv_path, query_fasta_path, receptor = chain2receptor[chain], CHMMAIRRa_bin = CHMMAIRRa_bin, uchime_mode = "balanced")
    global real_data_times    
    real_data_times = vcat(real_data_times, DataFrame("Method" => collect(keys(timing_tuple)),
    "Time (s)" => collect(values(timing_tuple)),
    "Dataset" => [dataset_name for i in 1:length(timing_tuple)]))
end
CSV.write(joinpath(benchmark_dir, "real_data_times.tsv"), real_data_times, delim = "\t")
times = vcat(sim_times, real_data_times)
# convert datasets to indices for easier plotting
times[!,"dataset_ind"] = [dataset2ind[el] for el in times[!,"Dataset"]];
CSV.write(joinpath(benchmark_dir, "times.tsv"), times, delim = "\t")