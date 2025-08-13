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
variation_method = "shazam"

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
IGH_shm_rates = [0.0, 0.05, 0.1, 0.2]

n_sequences = 10000
chimerism_rate = 0.05
chimeric_seqs_n, nonchimeric_seqs_n = Int(floor(n_sequences * chimerism_rate)), Int(floor(n_sequences * (1 - chimerism_rate)))
rng = MersenneTwister(random_seed)

ig_seqtype = "Ig"
V_seq = degap(reference_sets["OGRDB_IGHV_human_one_allele_per_gene"][2][1])
D_seq = degap(reference_sets["OGRDB_IGHD_human_one_allele_per_gene"][2][1])
J_seq = degap(reference_sets["OGRDB_IGHJ_human_one_allele_per_gene"][2][2])
test_sets = DataFrame()
for gene in ProgressBar(['d', 'j'])
    refset_name = "OGRDB_IGH$(uppercase(gene))_human_one_allele_per_gene"
    for ((shm1, shm2), location) in ProgressBar(Base.product(zip(IGH_shm_rates, IGH_shm_rates), ["random", "middle33"]))
        if (shm1 > 0.0) & ((gene == 'd') | (gene == 'j'))
            continue
        end
        @info "gene=$(gene), shm1=$(shm1), shm2=$(shm2), location=$(location)"
        global test_sets
        chimeric_names, chimeric_seqs, breakpoint_positions = random_chimeras(reference_sets[refset_name][1], degap.(reference_sets[refset_name][2]), rng,
            min_pos = name2positions[location][1],
            max_pos = name2positions[location][2],
            n = chimeric_seqs_n,
            variation_method = variation_method,
            min_shm1 = shm1,
            max_shm1 = shm1,
            min_shm2 = shm2,
            max_shm2 = shm2,
            gene = gene,
            V_seq = V_seq,
            D_seq = D_seq,
            J_seq = J_seq,
            ig_seqtype = ig_seqtype);
        nonchimeric_names, nonchimeric_seqs = random_nonchimeras(reference_sets[refset_name][1], degap.(reference_sets[refset_name][2]), rng,
            n = nonchimeric_seqs_n,
            variation_method = variation_method,
            min_shm = minimum([shm1, shm2]),
            max_shm = maximum([shm1, shm2]),
            gene = gene,
            V_seq = V_seq,
            D_seq = D_seq,
            J_seq = J_seq);
        chimeric_df = DataFrame(shm1 = shm1, shm2 = shm2, location = location, label = true, sequence_id = chimeric_names, sequence = chimeric_seqs, breakpoint_position = breakpoint_positions, refset_name = refset_name)
        nonchimeric_df = DataFrame(shm1 = shm1, shm2 = shm2, location = location, label = false, sequence_id = nonchimeric_names, sequence = nonchimeric_seqs, breakpoint_position = missing, refset_name = refset_name)
        test_sets = vcat(test_sets, chimeric_df, nonchimeric_df)
    end
end
CSV.write(joinpath(output_dir, "IGHDJ_shazam_test_sets_n=$(n_sequences)_chimerism_rate=$(chimerism_rate).tsv"), test_sets, delim = "\t")