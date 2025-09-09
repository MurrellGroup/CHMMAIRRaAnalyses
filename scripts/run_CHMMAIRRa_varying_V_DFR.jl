"""
Script to run CHMMAIRRa on all 5 real datasets in this paper.
"""

# load packages
using Pkg
Pkg.develop(path = "../../CHMMAIRRa.jl")

using Glob, CHMMAIRRa, DataFrames, CSV

min_V_DFRs = [0, 1, 2, 3, 4, 5]

timing_output_path = joinpath(@__DIR__, "../../outputs/timing_benchmarks/varying_V_DFR_times_with_internal_dedup.tsv")

# set up args
args = Dict()
args["subsample"] = 100_000
args["assignments"] = "final/filtered.tsv.gz"
args["V_fasta"] = "final/database/V.fasta"

igdiscover_dir = joinpath(@__DIR__, "../../data/igdiscover22/")

cd(igdiscover_dir)
TCR_D01_folders = filter(isdir, glob("GKH_TCR/*/D01/"))

IML3694_folders = ["IML369/IgM/IML3694_25cycle", "IML369/IgG/IML3694_post-vax"]

function julia_time_chmmairra(db_fasta_path::String, assignments_path::String, out_path::String; receptor::String = "TCR", detailed::Bool = false, min_V_DFR::Int = 0, subsample::Int = 10_000, disable_internal_dedup::Bool = false)
    t = @timed CHMMAIRRa.detect_chimeras_from_files(db_fasta_path, assignments_path, out_path,
            receptor = receptor, detailed = detailed, min_V_DFR = min_V_DFR, subsample = subsample, disable_internal_dedup = disable_internal_dedup)  
    return t.time
end

function precompile_chmmairra(db_fasta_path::String, assignments_path::String, out_path::String)
    println("Precompiling CHMMAIRRa")    
    julia_time_chmmairra(db_fasta_path, assignments_path, out_path, receptor = "IG", detailed = false, disable_internal_dedup = false)
    julia_time_chmmairra(db_fasta_path, assignments_path, out_path, receptor = "IG", detailed = true, disable_internal_dedup = false)
    julia_time_chmmairra(db_fasta_path, assignments_path, out_path, receptor = "TCR", detailed = false, disable_internal_dedup = false)
    julia_time_chmmairra(db_fasta_path, assignments_path, out_path, receptor = "TCR", detailed = true, disable_internal_dedup = false)
    println("Done precompiling CHMMAIRRa")
end

precompiled = false
detailed = false
times = DataFrame()
for (igdiscover_folder, min_V_DFR) in Base.product(TCR_D01_folders, min_V_DFRs)
    global times
    global precompiled

    cd(joinpath(igdiscover_dir, igdiscover_folder))
    args["out"] = "final/CHMMAIRRa_out.min-V-DFR=$(string(min_V_DFR)).internal-dedup=true.tsv.gz"
    if ! precompiled
        precompile_chmmairra(args["V_fasta"], args["assignments"], args["out"])
        precompiled = true
    end
    @info "Processing $igdiscover_folder with min_V_DFR = $min_V_DFR"
    t = julia_time_chmmairra(args["V_fasta"], args["assignments"], args["out"], receptor = "TCR", detailed = detailed, min_V_DFR = min_V_DFR, subsample = args["subsample"])
    out_df = CSV.read(args["out"], delim = "\t", DataFrame)
    chimeric_percent = sum(out_df.chimeric) / nrow(out_df) * 100
    times = vcat(times, DataFrame("case" => igdiscover_folder, "min_V_DFR" => min_V_DFR, "time" => t, "out" => args["out"], "detailed" => detailed, "subsample" => args["subsample"], "chimeric_percent" => chimeric_percent))
end

for (igdiscover_folder, min_V_DFR) in Base.product(IML3694_folders, min_V_DFRs)
    global times

    cd(joinpath(igdiscover_dir, igdiscover_folder))
    args["out"] = "final/CHMMAIRRa_out.min-V-DFR=$(string(min_V_DFR)).internal-dedup=true.tsv.gz"

    @info "Processing $igdiscover_folder with min_V_DFR = $min_V_DFR"
    t = julia_time_chmmairra(args["V_fasta"], args["assignments"], args["out"], receptor = "IG", detailed = detailed, min_V_DFR = min_V_DFR, subsample = args["subsample"])
    out_df = CSV.read(args["out"], delim = "\t", DataFrame)
    chimeric_percent = sum(out_df.chimeric) / nrow(out_df) * 100
    times = vcat(times, DataFrame("case" => igdiscover_folder, "min_V_DFR" => min_V_DFR, "time" => t, "out" => args["out"], "detailed" => detailed, "subsample" => args["subsample"], "chimeric_percent" => chimeric_percent))
end

CSV.write(timing_output_path, times, delim = "\t")