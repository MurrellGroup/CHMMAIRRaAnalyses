using Pkg
Pkg.activate("..")
# this makes sure we load the latest version of the code
using CHMMera, Random, CSV, DataFrames, Plots, PlotUtils, ProgressBars, StatsBase, JSON, Measures
include("../src/utils.jl")
include("../src/simulate.jl")

Threads.nthreads()

# INPUTS
# Only inputs needed are databases, which are provided with the repository
KITDB_dir = joinpath(@__DIR__, "../data/KI_TCR_DB_v0.0.1")  # See https://gkhlab.gitlab.io/tcr/sequences/
OGRDB_dir = joinpath(@__DIR__, "../data/OGRDB_human_IGH_9") # See https://ogrdb.airr-community.org/germline_sets/Homo%20sapiens
KI_TRB_dir = joinpath(KITDB_dir, "TRB")

# OUTPUTS
plots_dir = "../../outputs/plots/"


# mutation rates used for DB method in order to run on variable SHM rate simulated data
IG_DB_mutation_probabilities = [0.0, 0.005, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.25]
TCR_DB_mutation_probabilities = [0.005]
method2color = Dict("CHMMAIRRa BW" => "#BF40BF", "CHMMAIRRa DB" => "#0096FF", "USEARCH uchime2_ref" => "darkgreen", "VSEARCH uchime_ref" => "orange")
default_prior_probability = 0.05
random_seed = 888

all_methods_ROC_dir = joinpath(plots_dir, "ROCs/all_methods/")
CHMMAIRRa_only_ROC_dir = joinpath(plots_dir, "ROCs/CHMMAIRRa_only/")
rng = MersenneTwister(random_seed)
shm_method = "shazam"


#OGRDB_json = JSON.parsefile(joinpath(OGRDB_dir, "Homo_sapiens_IGH_VDJ_rev_9_ex.json"))
#OGRDB_refnames, OGRDB_refseqs = read_fasta(joinpath(OGRDB_dir, "V.fasta"));
## take only functional alleles
#OGRDB_functional_alleles = [el["label"] for el in OGRDB_json["GermlineSet"][1]["allele_descriptions"] if el["functional"] & occursin("IGHV", el["label"])]
#functional_inds = OGRDB_refnames .∈ Ref(OGRDB_functional_alleles)
#
#TRB_refnames, TRB_refseqs = read_fasta(joinpath(KI_TRB_dir, "V.fasta"));
#
## take one uniformally random allele per gene
#reference_sets = Dict("OGRDB_IGHV_human_one_allele_per_gene" => simulate_genotype(OGRDB_refnames[functional_inds], OGRDB_refseqs[functional_inds], rng))
#reference_sets["KI_TRBV_one_allele_per_gene"] = simulate_genotype(TRB_refnames, TRB_refseqs, rng)
#
## vary shm rates, choose either middlethird or random location for breakpoint
#IGH_shm_rates = [0.0, 0.05, 0.1, 0.2]
#test_sets = DataFrame()
#name2positions = Dict("random" => (0.0, 1.0), "middle90" => (0.05, .95), "middle80" => (0.1, .9), "middle60" => (0.2, 0.8), "middle33" => (1/3, 2/3))
#
#n_sequences = 10000
#chimerism_rate = 0.05
#chimeric_seqs_n, nonchimeric_seqs_n = Int(floor(n_sequences * chimerism_rate)), Int(floor(n_sequences * (1 - chimerism_rate)))
#rng = MersenneTwister(random_seed)
#
## IG simulations
#refset_name = "OGRDB_IGHV_human_one_allele_per_gene"
## for each combination of SHM rates and breakpoint location, generate chimeric and nonchimeric sequences
#for (shm1, shm2, location) in ProgressBar(Base.product(IGH_shm_rates, IGH_shm_rates, ["random", "middle33"]))
#    global test_sets
#    chimeric_names, chimeric_seqs, ref_names, ref_seqs, breakpoint_positions = random_chimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
#                                                                                            min_pos = name2positions[location][1],
#                                                                                            max_pos = name2positions[location][2],
#                                                                                            n = chimeric_seqs_n,
#                                                                                            min_shm1 = shm1,
#                                                                                            max_shm1 = shm1,
#                                                                                            min_shm2 = shm2,
#                                                                                            max_shm2 = shm2,
#                                                                                            shm_method = shm_method);
#    nonchimeric_names, nonchimeric_seqs = random_nonchimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
#                                                        n = nonchimeric_seqs_n,
#                                                        min_shm = maximum([shm1, shm2]),
#                                                        max_shm = maximum([shm1, shm2]),
#                                                        shm_method = shm_method);
#    chimeric_df = DataFrame(shm1 = shm1, shm2 = shm2, location = location, label = true, sequence_id = chimeric_names, sequence = chimeric_seqs, breakpoint_position = breakpoint_positions, refset_name = refset_name)
#    nonchimeric_df = DataFrame(shm1 = shm1, shm2 = shm2, location = location, label = false, sequence_id = nonchimeric_names, sequence = nonchimeric_seqs, breakpoint_position = missing, refset_name = refset_name)
#    test_sets = vcat(test_sets, chimeric_df, nonchimeric_df)
#end
#
## TRB simulations
#refset_name = "KI_TRBV_one_allele_per_gene"
#TRB_mutation_rates = [0.0, 0.001, 0.005, 0.02]
#for (TRB_mutation_rate, location) in ProgressBar(Base.product(TRB_mutation_rates, ["random", "middle33"]))
#    global test_sets
#    chimeric_names, chimeric_seqs, ref_names, ref_seqs, breakpoint_positions = random_chimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
#                                                                                            min_pos = name2positions[location][1],
#                                                                                            max_pos = name2positions[location][2],
#                                                                                            n = chimeric_seqs_n,
#                                                                                            min_shm1 = TRB_mutation_rate,
#                                                                                            max_shm1 = TRB_mutation_rate,
#                                                                                            min_shm2 = TRB_mutation_rate,
#                                                                                            max_shm2 = TRB_mutation_rate,
#                                                                                            shm_method = shm_method);
#    nonchimeric_names, nonchimeric_seqs = random_nonchimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
#                                                        n = nonchimeric_seqs_n,
#                                                        min_shm = TRB_mutation_rate,
#                                                        max_shm = TRB_mutation_rate,
#                                                        shm_method = shm_method);
#    chimeric_df = DataFrame(shm1 = TRB_mutation_rate, shm2 = TRB_mutation_rate, location = location, label = true, sequence_id = chimeric_names, sequence = chimeric_seqs, breakpoint_position = breakpoint_positions, refset_name = refset_name)
#    nonchimeric_df = DataFrame(shm1 = TRB_mutation_rate, shm2 = TRB_mutation_rate, location = location, label = false, sequence_id = nonchimeric_names, sequence = nonchimeric_seqs, breakpoint_position = missing, refset_name = refset_name)
#    test_sets = vcat(test_sets, chimeric_df, nonchimeric_df)
#end
#CSV.write("../../data/simulated/Vs_shazam_test_sets.tsv", test_sets, delim = "\t")


OGRDB_json = JSON.parsefile(joinpath(OGRDB_dir, "Homo_sapiens_IGH_VDJ_rev_9_ex.json"))
OGRDB_refnames, OGRDB_refseqs = read_fasta(joinpath(OGRDB_dir, "J.fasta"));
# take only functional alleles
OGRDB_functional_alleles = [el["label"] for el in OGRDB_json["GermlineSet"][1]["allele_descriptions"] if el["functional"] & occursin("IGHJ", el["label"])]
functional_inds = OGRDB_refnames .∈ Ref(OGRDB_functional_alleles)

TRB_refnames, TRB_refseqs = read_fasta(joinpath(KI_TRB_dir, "J.fasta"));

# take one uniformally random allele per gene
reference_sets = Dict("OGRDB_IGHJ_human_one_allele_per_gene" => simulate_genotype(OGRDB_refnames[functional_inds], OGRDB_refseqs[functional_inds], rng))
reference_sets["KI_TRBJ_one_allele_per_gene"] = simulate_genotype(TRB_refnames, TRB_refseqs, rng)


# vary shm rates, choose either middlethird or random location for breakpoint
IGH_shm_rates = [0.0, 0.05, 0.1, 0.2]
test_sets = DataFrame()
name2positions = Dict("random" => (0.0, 1.0), "middle90" => (0.05, .95), "middle80" => (0.1, .9), "middle60" => (0.2, 0.8), "middle33" => (1/3, 2/3))

n_sequences = 10000
chimerism_rate = 0.05

chimeric_seqs_n, nonchimeric_seqs_n = Int(floor(n_sequences * chimerism_rate)), Int(floor(n_sequences * (1 - chimerism_rate)))
rng = MersenneTwister(random_seed)

# IG simulations
refset_name = "OGRDB_IGHJ_human_one_allele_per_gene"
# for each combination of SHM rates and breakpoint location, generate chimeric and nonchimeric sequences
for (shm1, shm2, location) in ProgressBar(Base.product(IGH_shm_rates, IGH_shm_rates, ["random", "middle33"]))
    global test_sets
    chimeric_names, chimeric_seqs, ref_names, ref_seqs, breakpoint_positions = random_chimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
                                                                                            min_pos = name2positions[location][1],
                                                                                            max_pos = name2positions[location][2],
                                                                                            n = chimeric_seqs_n,
                                                                                            min_shm1 = shm1,
                                                                                            max_shm1 = shm1,
                                                                                            min_shm2 = shm2,
                                                                                            max_shm2 = shm2,
                                                                                            shm_method = shm_method);
    nonchimeric_names, nonchimeric_seqs = random_nonchimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
                                                        n = nonchimeric_seqs_n,
                                                        min_shm = maximum([shm1, shm2]),
                                                        max_shm = maximum([shm1, shm2]),
                                                        shm_method = shm_method);
    chimeric_df = DataFrame(shm1 = shm1, shm2 = shm2, location = location, label = true, sequence_id = chimeric_names, sequence = chimeric_seqs, breakpoint_position = breakpoint_positions, refset_name = refset_name)
    nonchimeric_df = DataFrame(shm1 = shm1, shm2 = shm2, location = location, label = false, sequence_id = nonchimeric_names, sequence = nonchimeric_seqs, breakpoint_position = missing, refset_name = refset_name)
    test_sets = vcat(test_sets, chimeric_df, nonchimeric_df)
end

# TRB simulations
refset_name = "KI_TRBJ_one_allele_per_gene"
TRB_mutation_rates = [0.0, 0.001, 0.005, 0.02]
for (TRB_mutation_rate, location) in ProgressBar(Base.product(TRB_mutation_rates, ["random", "middle33"]))
    global test_sets
    chimeric_names, chimeric_seqs, ref_names, ref_seqs, breakpoint_positions = random_chimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
                                                                                            min_pos = name2positions[location][1],
                                                                                            max_pos = name2positions[location][2],
                                                                                            n = chimeric_seqs_n,
                                                                                            min_shm1 = TRB_mutation_rate,
                                                                                            max_shm1 = TRB_mutation_rate,
                                                                                            min_shm2 = TRB_mutation_rate,
                                                                                            max_shm2 = TRB_mutation_rate,
                                                                                            shm_method = shm_method);
    nonchimeric_names, nonchimeric_seqs = random_nonchimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
                                                        n = nonchimeric_seqs_n,
                                                        min_shm = TRB_mutation_rate,
                                                        max_shm = TRB_mutation_rate,
                                                        shm_method = shm_method);
    chimeric_df = DataFrame(shm1 = TRB_mutation_rate, shm2 = TRB_mutation_rate, location = location, label = true, sequence_id = chimeric_names, sequence = chimeric_seqs, breakpoint_position = breakpoint_positions, refset_name = refset_name)
    nonchimeric_df = DataFrame(shm1 = TRB_mutation_rate, shm2 = TRB_mutation_rate, location = location, label = false, sequence_id = nonchimeric_names, sequence = nonchimeric_seqs, breakpoint_position = missing, refset_name = refset_name)
    test_sets = vcat(test_sets, chimeric_df, nonchimeric_df)
end
CSV.write("../../data/simulated/Js_shazam_test_sets.tsv", test_sets, delim = "\t")


OGRDB_json = JSON.parsefile(joinpath(OGRDB_dir, "Homo_sapiens_IGH_VDJ_rev_9_ex.json"))
OGRDB_refnames, OGRDB_refseqs = read_fasta(joinpath(OGRDB_dir, "D.fasta"));
# take only functional alleles
OGRDB_functional_alleles = [el["label"] for el in OGRDB_json["GermlineSet"][1]["allele_descriptions"] if el["functional"] & occursin("IGHD", el["label"])]
functional_inds = OGRDB_refnames .∈ Ref(OGRDB_functional_alleles)

TRB_refnames, TRB_refseqs = read_fasta(joinpath(KI_TRB_dir, "D.fasta"));

# take one uniformally random allele per gene
reference_sets = Dict("OGRDB_IGHD_human_one_allele_per_gene" => simulate_genotype(OGRDB_refnames[functional_inds], OGRDB_refseqs[functional_inds], rng))
reference_sets["KI_TRBD_one_allele_per_gene"] = simulate_genotype(TRB_refnames, TRB_refseqs, rng)


# vary shm rates, choose either middlethird or random location for breakpoint
IGH_shm_rates = [0.0, 0.05, 0.1, 0.2]
test_sets = DataFrame()
name2positions = Dict("random" => (0.0, 1.0), "middle90" => (0.05, .95), "middle80" => (0.1, .9), "middle60" => (0.2, 0.8), "middle33" => (1/3, 2/3))

n_sequences = 10000
chimerism_rate = 0.05

chimeric_seqs_n, nonchimeric_seqs_n = Int(floor(n_sequences * chimerism_rate)), Int(floor(n_sequences * (1 - chimerism_rate)))
rng = MersenneTwister(random_seed)

# IG simulations
refset_name = "OGRDB_IGHD_human_one_allele_per_gene"
# for each combination of SHM rates and breakpoint location, generate chimeric and nonchimeric sequences
for (shm1, shm2, location) in ProgressBar(Base.product(IGH_shm_rates, IGH_shm_rates, ["random", "middle33"]))
    global test_sets
    chimeric_names, chimeric_seqs, ref_names, ref_seqs, breakpoint_positions = random_chimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
                                                                                            min_pos = name2positions[location][1],
                                                                                            max_pos = name2positions[location][2],
                                                                                            n = chimeric_seqs_n,
                                                                                            min_shm1 = shm1,
                                                                                            max_shm1 = shm1,
                                                                                            min_shm2 = shm2,
                                                                                            max_shm2 = shm2,
                                                                                            shm_method = shm_method);
    nonchimeric_names, nonchimeric_seqs = random_nonchimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
                                                        n = nonchimeric_seqs_n,
                                                        min_shm = maximum([shm1, shm2]),
                                                        max_shm = maximum([shm1, shm2]),
                                                        shm_method = shm_method);
    chimeric_df = DataFrame(shm1 = shm1, shm2 = shm2, location = location, label = true, sequence_id = chimeric_names, sequence = chimeric_seqs, breakpoint_position = breakpoint_positions, refset_name = refset_name)
    nonchimeric_df = DataFrame(shm1 = shm1, shm2 = shm2, location = location, label = false, sequence_id = nonchimeric_names, sequence = nonchimeric_seqs, breakpoint_position = missing, refset_name = refset_name)
    test_sets = vcat(test_sets, chimeric_df, nonchimeric_df)
end

# TRB simulations
refset_name = "KI_TRBD_one_allele_per_gene"
TRB_mutation_rates = [0.0, 0.001, 0.005, 0.02]
for (TRB_mutation_rate, location) in ProgressBar(Base.product(TRB_mutation_rates, ["random", "middle33"]))
    global test_sets
    chimeric_names, chimeric_seqs, ref_names, ref_seqs, breakpoint_positions = random_chimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
                                                                                            min_pos = name2positions[location][1],
                                                                                            max_pos = name2positions[location][2],
                                                                                            n = chimeric_seqs_n,
                                                                                            min_shm1 = TRB_mutation_rate,
                                                                                            max_shm1 = TRB_mutation_rate,
                                                                                            min_shm2 = TRB_mutation_rate,
                                                                                            max_shm2 = TRB_mutation_rate,
                                                                                            shm_method = shm_method);
    nonchimeric_names, nonchimeric_seqs = random_nonchimeras(reference_sets[refset_name][1], reference_sets[refset_name][2], rng,
                                                        n = nonchimeric_seqs_n,
                                                        min_shm = TRB_mutation_rate,
                                                        max_shm = TRB_mutation_rate,
                                                        shm_method = shm_method);
    chimeric_df = DataFrame(shm1 = TRB_mutation_rate, shm2 = TRB_mutation_rate, location = location, label = true, sequence_id = chimeric_names, sequence = chimeric_seqs, breakpoint_position = breakpoint_positions, refset_name = refset_name)
    nonchimeric_df = DataFrame(shm1 = TRB_mutation_rate, shm2 = TRB_mutation_rate, location = location, label = false, sequence_id = nonchimeric_names, sequence = nonchimeric_seqs, breakpoint_position = missing, refset_name = refset_name)
    test_sets = vcat(test_sets, chimeric_df, nonchimeric_df)
end
CSV.write("../../data/simulated/Ds_shazam_test_sets.tsv", test_sets, delim = "\t")
