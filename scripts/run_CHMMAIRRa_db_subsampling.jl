#!/bin/bash

################################################
#
# Run IgBLAST and CHMMAIRRa commands to generate data for the subsampling test in figure 4
#
################################################

# run 20 different seeds per database size
# D01 D04 EUR, D07 SAS, D14 EAS, D11 AFR

# ran with IgDiscover22 v1.0.2

using Pkg
Pkg.activate(".")
Pkg.develop(path = "../CHMMAIRRa.jl")
using CHMMAIRRa, Glob, CSV, DataFrames, FASTX, Random, StatsBase, CodecZlib

include("./src/utils.jl")

threads = 16
seeds = [1, 2, 3]
chain2sizes = Dict("IGH" => [30, 40],
                   "TRA" => [30, 40],
                   "TRB" => [40, 50],
                   "TRG" => [4, 6])

# set up paths to databases

IMGT_VQUEST_dir = "/home/mchernys/references/IMGT/V-QUEST/homo_sapiens/202430-2/"
human_gl_aux = "$(IMGT_VQUEST_dir)/human_gl.aux"
IMGT_dbs = Dict("IGH" => "$(IMGT_VQUEST_dir)/IG/F+P/H",
                "TRA" => "$(IMGT_VQUEST_dir)/TR/F+P/A",
                "TRB" => "$(IMGT_VQUEST_dir)/TR/F+P/B",
                "TRG" => "$(IMGT_VQUEST_dir)/TR/F+P/G")

IMGT_F_dbs = Dict("IGH" => "$(IMGT_VQUEST_dir)/IG/F/H",
                  "TRA" => "$(IMGT_VQUEST_dir)/TR/F/A",
                  "TRB" => "$(IMGT_VQUEST_dir)/TR/F/B",
                  "TRG" => "$(IMGT_VQUEST_dir)/TR/F/G")

KITDB_dir = "/home/mchernys/references/KITDB"
KI_TCR_dbs = Dict("TRA" => "$(KITDB_dir)/TRA",
                  "TRB" => "$(KITDB_dir)/TRB",
                  "TRG" => "$(KITDB_dir)/TRG")

OGRDB_dir = "/home/mchernys/references/OGRDB/human/IGH/9/"

# set up paths to IgDiscover runs

igdiscover_dir = "/home/mchernys/ben/chimera_detection/data/igdiscover22/"
igdiscover_run_dirs = Dict("IGH" => ["$(igdiscover_dir)IML3694/runs/IgG/IML3694_acute",
                                     "$(igdiscover_dir)IML3694/runs/IgG/IML3694_convalescent",
                                     "$(igdiscover_dir)IML3694/runs/IgG/IML3694_post-vax",
                                     "$(igdiscover_dir)IML3694/runs/IgM/IML3694_25cycle",
                                     "$(igdiscover_dir)IML3695/runs/IgG/IML3695_acute",
                                     "$(igdiscover_dir)IML3695/runs/IgG/IML3695_convalescent",
                                     "$(igdiscover_dir)IML3695/runs/IgG/IML3695_post-vax",
                                     "$(igdiscover_dir)IML3695/runs/IgM/IML3695_25cycle"],

                            "TRA" => ["$(igdiscover_dir)GKH_TCR/TRA/D01",
                                      "$(igdiscover_dir)GKH_TCR/TRA/D04",
                                      "$(igdiscover_dir)GKH_TCR/TRA/D07",
                                      "$(igdiscover_dir)GKH_TCR/TRA/D14",
                                      "$(igdiscover_dir)GKH_TCR/TRA/D11"],

                            "TRB" => ["$(igdiscover_dir)GKH_TCR/TRB/D01",
                                      "$(igdiscover_dir)GKH_TCR/TRB/D04",
                                      "$(igdiscover_dir)GKH_TCR/TRB/D07",
                                      "$(igdiscover_dir)GKH_TCR/TRB/D14",
                                      "$(igdiscover_dir)GKH_TCR/TRB/D11"],

                            "TRG" => ["$(igdiscover_dir)GKH_TCR/TRG/D01",
                                      "$(igdiscover_dir)GKH_TCR/TRG/D04",
                                      "$(igdiscover_dir)GKH_TCR/TRG/D07",
                                      "$(igdiscover_dir)GKH_TCR/TRG/D14",
                                      "$(igdiscover_dir)GKH_TCR/TRG/D11"])

# set up CHMMAIRRa args
args = Dict{String, Any}("align-database" => true,
            "detailed" => true,
            "count-chimeric-segments" => true)

function create_subsampled_db(igdiscover_dir::String, db_dir::String, db_size::Int, seed::Int)
    cd(igdiscover_dir)
    refnames, refseqs = read_fasta("final/database/V.fasta")
    Random.seed!(seed)
    inds = sample(collect(1:length(refseqs)), db_size, replace=false)
    if ! isdir(db_dir)
        mkdir(db_dir)
    end
    write_fasta(joinpath(db_dir, "V.fasta"), refseqs[inds], seq_names = refnames[inds])
    cp("final/database/D.fasta", joinpath(db_dir, "D.fasta"), force = true)
    cp("final/database/J.fasta", joinpath(db_dir, "J.fasta"), force = true)
    cp("database/human_gl.aux", joinpath(db_dir, "human_gl.aux"), force = true)
    return db_dir
end


for chain in ["IGH", "TRA", "TRB", "TRG"]
    @info "Processing $(chain)"
    args["receptor"] = chain == "IGH" ? "IG" : "TCR"
    # PG subsets
    for igdiscover_run_dir in igdiscover_run_dirs[chain]
        cd(igdiscover_run_dir)
        for size in chain2sizes[chain]
            for seed in seeds
                @info "Processing $(igdiscover_run_dir) PG.$(size)Vs.$(seed)seed"
                PG_subset_dir = "$(size)Vs-$(seed)seed-database"
                #create_subsampled_db(igdiscover_run_dir, PG_subset_dir, size, seed)
                #run(`conda run -n igdiscover igdiscover igblastwrap --sequence-type $(sequence_type) --aux $(human_gl_aux)--threads $(threads) $(PG_subset_dir) reads/sequences_50000subsample.fasta.gz | pigz -c > final/$(size)Vs-$(seed)seed-airr.tsv.gz`)
                #run(`conda run -n igdiscover igdiscover augment --read-cdr3 $(PG_subset_dir) final/$(size)Vs-$(seed)seed-airr.tsv.gz | pigz -c > final/$(size)Vs-$(seed)seed-assigned.tsv.gz`)
                #run(`conda run -n igdiscover igdiscover filter final/$(size)Vs-$(seed)seed-assigned.tsv.gz | pigz -c > final/$(size)Vs-$(seed)seed-filtered.tsv.gz`)
                args["V_fasta"] = "$(PG_subset_dir)/V.fasta"
                args["assignments"] = "final/$(size)Vs-$(seed)seed-filtered.tsv.gz"
                args["out"] = "final/CHMMAIRRa_out.PG.$(size)Vs.$(seed)seed.tsv.gz"
                args["chimeric-alignments"] = "final/CHMMAIRRa_out.PG.$(size)Vs.$(seed)seed.fasta"
                CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = args["receptor"],
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"])
            end
        end

        # PG
        @info "Processing $(igdiscover_run_dir) PG"
        #run(`conda run -n igdiscover igdiscover igblastwrap --sequence-type $(sequence_type) --aux $(human_gl_aux) --threads $(threads) final/database/ reads/sequences_50000subsample.fasta.gz | pigz -c > final/PG-airr.tsv.gz;`)
        #run(`conda run -n igdiscover igdiscover augment --read-cdr3 final/database/ final/PG-airr.tsv.gz | pigz -c > final/PG-assigned.tsv.gz;`)
        #run(`conda run -n igdiscover igdiscover filter final/PG-assigned.tsv.gz | pigz -c > final/PG-filtered.tsv.gz;`)
        args["V_fasta"] = "final/database/V.fasta"
        args["assignments"] = "final/PG-filtered.tsv.gz"
        args["out"] = "final/CHMMAIRRa_out.PG.tsv.gz"
        args["chimeric-alignments"] = "final/CHMMAIRRa_out.PG.fasta"
        CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                    receptor = args["receptor"],
                    align_database = args["align-database"],
                    detailed = args["detailed"],
                    count_chimeric_segments = args["count-chimeric-segments"],
                    chimeric_alignments = args["chimeric-alignments"])
        ## IMGT db comparison
        #@info "Processing $(igdiscover_run_dir) IMGT db"
        ##run(`conda run -n igdiscover igdiscover igblastwrap --sequence-type $(sequence_type) --aux $(human_gl_aux) --threads $(threads) $(IMGT_dbs[chain]) reads/sequences_50000subsample.fasta.gz | pigz -c > final/IMGT-V-QUEST-F+P-airr.tsv.gz;`)
        ##run(`conda run -n igdiscover igdiscover augment --read-cdr3 $(IMGT_dbs[chain]) final/IMGT-V-QUEST-F+P-airr.tsv.gz | pigz -c > final/IMGT-V-QUEST-F+P-assigned.tsv.gz;`)
        ##run(`conda run -n igdiscover igdiscover filter final/IMGT-V-QUEST-F+P-assigned.tsv.gz | pigz -c > final/IMGT-V-QUEST-F+P-filtered.tsv.gz;`)
        #args["V_fasta"] = "$(IMGT_dbs[chain])/V.fasta"
        #args["assignments"] = "final/IMGT-V-QUEST-F+P-filtered.tsv.gz"
        #args["out"] = "final/CHMMAIRRa_out.IMGT-V-QUEST-F+P.tsv.gz"
        #args["chimeric-alignments"] = "final/CHMMAIRRa_out.IMGT-V-QUEST-F+P.fasta"
        #CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
        #            receptor = args["receptor"],
        #            align_database = args["align-database"],
        #            detailed = args["detailed"],
        #            count_chimeric_segments = args["count-chimeric-segments"],
        #            chimeric_alignments = args["chimeric-alignments"])

        # IMGT F db comparison
        @info "Processing $(igdiscover_run_dir) IMGT F db"
        #run(`conda run -n igdiscover igdiscover igblastwrap --sequence-type $(sequence_type) --aux $(human_gl_aux) --threads $(threads) $(IMGT_F_dbs[chain]) reads/sequences_50000subsample.fasta.gz | pigz -c > final/IMGT-V-QUEST-F-airr.tsv.gz;`)
        #run(`conda run -n igdiscover igdiscover augment --read-cdr3 $(IMGT_F_dbs[chain]) final/IMGT-V-QUEST-F-airr.tsv.gz | pigz -c > final/IMGT-V-QUEST-F-assigned.tsv.gz;`)
        #run(`conda run -n igdiscover igdiscover filter final/IMGT-V-QUEST-F-assigned.tsv.gz | pigz -c > final/IMGT-V-QUEST-F-filtered.tsv.gz;`)
        args["V_fasta"] = "$(IMGT_dbs[chain])/V.fasta"
        args["assignments"] = "final/IMGT-V-QUEST-F-filtered.tsv.gz"
        args["out"] = "final/CHMMAIRRa_out.IMGT-V-QUEST-F.tsv.gz"
        args["chimeric-alignments"] = "final/CHMMAIRRa_out.IMGT-V-QUEST-F.fasta"
        CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                    receptor = args["receptor"],
                    align_database = args["align-database"],
                    detailed = args["detailed"],
                    count_chimeric_segments = args["count-chimeric-segments"],
                    chimeric_alignments = args["chimeric-alignments"])

        if chain == "IGH"
            @info "Processing $(igdiscover_run_dir) OGRDB"
            #run_igblastwrap_from_files(OGRDB_dir, "reads/sequences_50000subsample.fasta.gz", "final/OGRDB-airr.tsv.gz", sequence_type = "Ig", human_gl_aux = human_gl_aux, threads = threads)
            #run_augment_from_files(OGRDB_dir, "final/OGRDB-airr.tsv.gz", "final/OGRDB-assigned.tsv.gz")
            #run_filter_from_files("final/OGRDB-assigned.tsv.gz", "final/OGRDB-filtered.tsv.gz")
            args["V_fasta"] = "$(OGRDB_dir)/V.fasta"
            args["assignments"] = "final/OGRDB-filtered.tsv.gz"
            args["out"] = "final/CHMMAIRRa_out.OGRDB.tsv.gz"
            args["chimeric-alignments"] = "final/CHMMAIRRa_out.OGRDB.fasta"
            CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = args["receptor"],
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"])
        end

        if chain != "IGH"
            # for KITDB
            @info "Processing $(igdiscover_run_dir) KITDB"
            #run(`conda run -n igdiscover igdiscover igblastwrap --sequence-type $(sequence_type) --aux $(human_gl_aux) --threads $(threads) $(KI_TCR_dbs[chain]) reads/sequences_50000subsample.fasta.gz | pigz -c > final/KITDB-airr.tsv.gz;`)
            #run(`conda run -n igdiscover igdiscover augment --read-cdr3 $(KI_TCR_dbs[chain]) final/KITDB-airr.tsv.gz | pigz -c > final/KITDB-assigned.tsv.gz;`)
            #run(`conda run -n igdiscover igdiscover filter final/KITDB-assigned.tsv.gz | pigz -c > final/KITDB-filtered.tsv.gz;`)
            args["V_fasta"] = "$(KI_TCR_dbs[chain])/V.fasta"
            args["assignments"] = "final/KITDB-filtered.tsv.gz"
            args["out"] = "final/CHMMAIRRa_out.KITDB.tsv.gz"
            args["chimeric-alignments"] = "final/CHMMAIRRa_out.KITDB.fasta"
            CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = args["receptor"],
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"])
        end
    end
end

# set up paths to databases

IMGT_VQUEST_dir = "/home/mchernys/references/IMGT/V-QUEST/homo_sapiens/202430-2/"
IMGT_dbs = Dict("IGH" => "$(IMGT_VQUEST_dir)/IG/F+P/H",
                "TRA" => "$(IMGT_VQUEST_dir)/TR/F+P/A",
                "TRB" => "$(IMGT_VQUEST_dir)/TR/F+P/B",
                "TRG" => "$(IMGT_VQUEST_dir)/TR/F+P/G")

IMGT_F_dbs = Dict("IGH" => "$(IMGT_VQUEST_dir)/IG/F/H",
                  "TRA" => "$(IMGT_VQUEST_dir)/TR/F/A",
                  "TRB" => "$(IMGT_VQUEST_dir)/TR/F/B",
                  "TRG" => "$(IMGT_VQUEST_dir)/TR/F/G")

# set up paths to IgDiscover runs

igdiscover_dir = "/home/mchernys/ben/chimera_detection/data/igdiscover22/"

# set up CHMMAIRRa args
args = Dict{String, Any}("align-database" => true,
            "detailed" => true,
            "count-chimeric-segments" => true)

# get recombfreqplot for D04 for supplementary figure 1
cd("$(igdiscover_dir)GKH_TCR/TRG/D04")
# IMGT
chain = "TRG"
args["V_fasta"] = "$(IMGT_dbs[chain])/V.fasta"
args["assignments"] = "final/IMGT-V-QUEST-F-filtered.tsv.gz"
args["out"] = "final/CHMMAIRRa_out.IMGT-V-QUEST-F.tsv.gz"
args["chimeric-alignments"] = "final/CHMMAIRRa_out.IMGT-V-QUEST-F.fasta"
CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
            receptor = "TCR",
            align_database = args["align-database"],
            detailed = args["detailed"],
            count_chimeric_segments = args["count-chimeric-segments"],
            chimeric_alignments = args["chimeric-alignments"],
            recombfreqplot = "final/TRG_D04.recombfreqplot.IMGT-VQUEST-F.svg")

# PG
args["V_fasta"] = "final/database/V.fasta"
args["assignments"] = "final/PG-filtered.tsv.gz"
args["out"] = "final/CHMMAIRRa_out.PG.tsv.gz"
args["chimeric-alignments"] = "final/CHMMAIRRa_out.PG.fasta"
CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
            receptor = "TCR",
            align_database = args["align-database"],
            detailed = args["detailed"],
            count_chimeric_segments = args["count-chimeric-segments"],
            chimeric_alignments = args["chimeric-alignments"],
            recombfreqplot = "final/TRG_D04.recombfreqplot.PG.svg")