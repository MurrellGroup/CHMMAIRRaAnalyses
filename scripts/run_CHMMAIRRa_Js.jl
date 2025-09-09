"""
Script to run CHMMAIRRa on all 5 real datasets in this paper.
"""

# load packages
using Pkg
Pkg.develop(path = "../../CHMMAIRRa.jl")

using Glob, CHMMAIRRa, DataFrames, CSV

# set up args
args = Dict()
args["subsample"] = 10_000
args["chimeric-alignments"] = "final/J_chimeric_alignments.fasta"
args["assignments"] = "final/filtered.tsv.gz"
args["V_fasta"] = "final/database/V.fasta"
args["J_fasta"] = "final/database/J.fasta"
args["recombfreqplot"] = "final/J_recombfreqplot.tsv"
args["out"] = "final/CHMMAIRRa_out.J.tsv.gz"
igdiscover_dir = joinpath(@__DIR__, "../../data/igdiscover22/")

cd(igdiscover_dir)
TCR_D01_folders = ["GKH_TCR/TRA/D01/", "GKH_TCR/TRB/D01/", "GKH_TCR/TRG/D01/"]

IML3694_folders = ["IML369/IgM/IML3694_25cycle/", "IML369/IgG/IML3694_post-vax/"]

for igdiscover_folder in TCR_D01_folders

    cd(joinpath(igdiscover_dir, igdiscover_folder))
    @info "Processing $igdiscover_folder"

    CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                                        J_fasta = args["J_fasta"],
                                        receptor = "TCR",
                                        detailed = true, 
                                        chimeric_alignments = args["chimeric-alignments"],
                                        recombfreqplot = args["recombfreqplot"],
                                        subsample = args["subsample"])
end

for igdiscover_folder in IML3694_folders

    cd(joinpath(igdiscover_dir, igdiscover_folder))
    
    CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                                        J_fasta = args["J_fasta"],
                                        receptor = "IG",
                                        detailed = true, 
                                        chimeric_alignments = args["chimeric-alignments"],
                                        recombfreqplot = args["recombfreqplot"],
                                        subsample = args["subsample"])
end