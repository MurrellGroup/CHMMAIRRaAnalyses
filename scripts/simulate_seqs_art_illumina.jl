cd("/home/mchernys/ben/chimera_detection/CHMMAIRRaAnalyses/scripts")
using Pkg
Pkg.activate("..")

# this makes sure we load the latest version of the code
using Random, CSV, DataFrames, Plots, PlotUtils, ProgressBars, StatsBase, JSON, Measures, BioAlignments, IgBLAST
include("../src/utils.jl")
include("../src/simulate.jl")

# INPUTS
# Only inputs needed are databases, which are provided with the repository
KITDB_dir = joinpath(@__DIR__, "../data/KI_TCR_DB_v0.0.1")  # See https://gkhlab.gitlab.io/tcr/sequences/
OGRDB_dir = joinpath(@__DIR__, "../data/OGRDB_human_IGH_9") # See https://ogrdb.airr-community.org/germline_sets/Homo%20sapiens
KI_TRB_dir = joinpath(KITDB_dir, "TRB")

# OUTPUTS
output_dir = "../../data/simulated/"

random_seed = 888




rng = MersenneTwister(random_seed)
variation_method = "art_illumina"

OGRDB_json = JSON.parsefile(joinpath(OGRDB_dir, "Homo_sapiens_IGH_VDJ_rev_9_ex.json"))

# generate synthetic IGHV, IGHD, IGHJ, TRBV, TRBD, TRBJ genotypes from OGRDB and KI
# one allele per gene
reference_sets = Dict()
for gene in ["V", "D", "J"]
    OGRDB_refnames, OGRDB_refseqs = read_fasta(joinpath(OGRDB_dir, "$(gene).fasta"));
    # take only functional alleles
    OGRDB_functional_alleles = [el["label"] for el in OGRDB_json["GermlineSet"][1]["allele_descriptions"] if el["functional"] & occursin("IGH$(gene)", el["label"])]
    functional_inds = OGRDB_refnames .∈ Ref(OGRDB_functional_alleles)
    TRB_refnames, TRB_refseqs = read_fasta(joinpath(KI_TRB_dir, "$(gene).fasta"));
    reference_sets["OGRDB_IGH$(gene)_human_one_allele_per_gene"] = simulate_genotype(OGRDB_refnames[functional_inds], OGRDB_refseqs[functional_inds], rng)
    reference_sets["KI_TRB$(gene)_one_allele_per_gene"] = simulate_genotype(TRB_refnames, TRB_refseqs, rng)
    write_fasta(joinpath(output_dir, "IGH$(gene)_one_allele_per_gene.fasta"), reference_sets["OGRDB_IGH$(gene)_human_one_allele_per_gene"][2], seq_names = reference_sets["OGRDB_IGH$(gene)_human_one_allele_per_gene"][1])
    write_fasta(joinpath(output_dir, "TRB$(gene)_one_allele_per_gene.fasta"), reference_sets["KI_TRB$(gene)_one_allele_per_gene"][2], seq_names = reference_sets["KI_TRB$(gene)_one_allele_per_gene"][1])
end

name2positions = Dict("random" => (0.0, 1.0), "middle90" => (0.05, .95), "middle80" => (0.1, .9), "middle60" => (0.2, 0.8), "middle33" => (1/3, 2/3))

n_sequences = 1000
chimerism_rate = 0.05
chimeric_seqs_n, nonchimeric_seqs_n = Int(floor(n_sequences * chimerism_rate)), Int(floor(n_sequences * (1 - chimerism_rate)))
rng = MersenneTwister(random_seed)

include("../src/utils.jl")
include("../src/simulate.jl")
V_seq = degap(reference_sets["KI_TRBV_one_allele_per_gene"][2][1])
D_seq = degap(reference_sets["KI_TRBD_one_allele_per_gene"][2][1])
J_seq = degap(reference_sets["KI_TRBJ_one_allele_per_gene"][2][2])
test_sets = DataFrame()
for gene in ProgressBar(["v", "d", "j"])
    @info "Simulating $(gene) sequences"
    refset_name = "KI_TRB$(uppercase(gene))_one_allele_per_gene"
    for location in ["random", "middle33"]
        @info "location=$(location)"
        global test_sets
        chimeric_names, chimeric_seqs, breakpoint_positions = random_chimeras(reference_sets[refset_name][1], degap.(reference_sets[refset_name][2]), rng,
            min_pos = name2positions[location][1],
            max_pos = name2positions[location][2],
            n = chimeric_seqs_n,
            variation_method = variation_method,
            gene = gene,
            V_seq = V_seq,
            D_seq = D_seq,
            J_seq = J_seq);
        nonchimeric_names, nonchimeric_seqs = random_nonchimeras(reference_sets[refset_name][1], degap.(reference_sets[refset_name][2]), rng,
            n = nonchimeric_seqs_n,
            variation_method = variation_method,
            gene = gene,
            V_seq = V_seq,
            D_seq = D_seq,
            J_seq = J_seq);
        chimeric_df = DataFrame(location = location, label = true, sequence_id = chimeric_names, sequence = chimeric_seqs, breakpoint_position = breakpoint_positions, refset_name = refset_name)
        nonchimeric_df = DataFrame(location = location, label = false, sequence_id = nonchimeric_names, sequence = nonchimeric_seqs, breakpoint_position = missing, refset_name = refset_name)
        test_sets = vcat(test_sets, chimeric_df, nonchimeric_df)
    end
end
CSV.write(joinpath(output_dir, "TRB_art_illumina_test_sets.tsv"), test_sets, delim = "\t")