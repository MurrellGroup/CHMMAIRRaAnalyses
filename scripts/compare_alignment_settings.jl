"""
Script to run CHMMAIRRa 5 datasets across all adaptive immune receptor chains using 3 database alignment methods.
Alignment methods:
- mafft_default: default mafft settings
- mafft_localpair: mafft with localpair option
- muscle_default: default muscle settings
"""

# load packages
using Pkg
Pkg.activate("..")
include("../src/utils.jl")

using CHMMAIRRa, Glob, CSV, DataFrames


# set up args
args = Dict()
args["align-database"] = false
args["detailed"] = true
args["count-chimeric-segments"] = true
args["subsample"] = 100_000
args["assignments"] = "final/filtered.tsv.gz"

igdiscover_dir = joinpath(@__DIR__, "../../data/igdiscover22/")

cd(igdiscover_dir)
TCR_D01_folders = filter(isdir, glob("GKH_TCR/*/D01/"))

IML3694_folders = ["IML369/IgM/IML3694_25cycle", "IML369/IgG/IML3694_post-vax"]

function align_database(sequence_ids::Vector{String}, sequences::Vector{String}, alignment_setting::String)
    if alignment_setting == "mafft_default"
        return mafft_wrapper(sequences, sequence_ids)
    elseif alignment_setting == "mafft_localpair"
        return mafft_wrapper(sequences, sequence_ids, parameters = ["--localpair", "--maxiterate", "1000"])
    elseif alignment_setting == "muscle_default"
        return muscle_wrapper(sequences, sequence_ids)
    end
end


for (igdiscover_folder, alignment_setting) in Base.product(TCR_D01_folders, ["mafft_localpair", "muscle_default"])
    cd(joinpath(igdiscover_dir, igdiscover_folder))

    args["V_fasta"] = "final/database/$(alignment_setting)_V.fasta"
    args["chimeric-alignments"] = "final/chimeric_$(alignment_setting).fasta"
    args["out"] = "final/CHMMAIRRa_out.$(alignment_setting).tsv.gz"
    V_names, V_seqs = read_fasta("final/database/V.fasta")
    aligned_V_names, aligned_V_seqs = align_database(V_names, V_seqs, alignment_setting)
    write_fasta(args["V_fasta"], aligned_V_seqs, seq_names = aligned_V_names)
    
    @info "Processing $igdiscover_folder with $alignment_setting"
    @time CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = "TCR",
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"])
end

for (igdiscover_folder, alignment_setting) in Base.product(IML3694_folders, ["mafft_localpair", "muscle_default"])
    cd(joinpath(igdiscover_dir, igdiscover_folder))

    args["V_fasta"] = "final/database/$(alignment_setting)_V.fasta"
    args["chimeric-alignments"] = "final/chimeric_$(alignment_setting).fasta"
    args["out"] = "final/CHMMAIRRa_out.$(alignment_setting).tsv.gz"
    V_names, V_seqs = read_fasta("final/database/V.fasta")
    aligned_V_names, aligned_V_seqs = align_database(V_names, V_seqs, alignment_setting)
    write_fasta(args["V_fasta"], aligned_V_seqs, seq_names = aligned_V_names)
    
    @info "Processing $igdiscover_folder with $alignment_setting"
    @time CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = "IG",
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"],
                        subsample = args["subsample"],
                        )
end