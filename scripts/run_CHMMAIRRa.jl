"""
Script to run CHMMAIRRa on all 5 real datasets in this paper.
"""

# load packages
using Pkg
Pkg.activate("..")
using CHMMAIRRa, Glob, CSV, DataFrames


# set up args
args = Dict()
args["align-database"] = true
args["detailed"] = true
args["count-chimeric-segments"] = true
args["chimeric-airr"] = "final/chimeric_22-1-25.tsv.gz"
args["chimeric-alignments"] = "final/chimeric_22-1-25.fasta"
args["V_fasta"] = "final/database/V.fasta"
args["assignments"] = "final/filtered.tsv.gz"
args["out"] = "final/CHMMAIRRa_out_22-1-25.tsv.gz"

igdiscover_dir = "../../data/igdiscover22/"

cd(igdiscover_dir)

PCR_conditions_igdiscover_folders = filter(isdir, glob("PCR_conditions/D*"))

PRJNA308641_igdiscover_folders = filter(isdir, glob("PRJNA308641/Ig[MG]/SRR*"))

IML369_igdiscover_folders = filter(isdir, glob("IML369/Ig[MG]/IML*"))

PRJNA300878_igdiscover_folders = filter(isdir, glob("PRJNA300878/TR[AB]/*"))

GKH_TCR_igdiscover_folders = vcat(filter(isdir, glob("GKH_TCR/TRA/D*")),
                                sort(filter(isdir, glob("GKH_TCR/TRB/D*"))),
                                sort(filter(isdir, glob("GKH_TCR/TRD/D*"))),
                                sort(filter(isdir, glob("GKH_TCR/TRG/D*"))))

for PCR_conditions_igdiscover_folder in PCR_conditions_igdiscover_folders
    cd(joinpath(igdiscover_dir, PCR_conditions_igdiscover_folder))
    @info "Processing $PCR_conditions_igdiscover_folder"
    @time CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = "IG",
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"],
                        )
end

for PRJNA308641_igdiscover_folder in PRJNA308641_igdiscover_folders
    cd(joinpath(igdiscover_dir, PRJNA308641_igdiscover_folder))
    @info "Processing $PRJNA308641_igdiscover_folder"
    @time CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = "IG",
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"],
                        )
end

for IML369_igdiscover_folder in IML369_igdiscover_folders
    cd(joinpath(igdiscover_dir, IML369_igdiscover_folder))
    @info "Processing $IML369_igdiscover_folder"
    @time CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = "IG",
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"],
                        )
end


for PRJNA300878_igdiscover_folder in PRJNA300878_igdiscover_folders
    cd(joinpath(igdiscover_dir, PRJNA300878_igdiscover_folder))
    @info "Processing $PRJNA300878_igdiscover_folder"
    @time CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = "TCR",
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"],
                        )
end

for GKH_TCR_igdiscover_folder in GKH_TCR_igdiscover_folders
    cd(joinpath(igdiscover_dir, GKH_TCR_igdiscover_folder))
    @info "Processing $GKH_TCR_igdiscover_folder"
    @time CHMMAIRRa.detect_chimeras_from_files(args["V_fasta"], args["assignments"], args["out"],
                        receptor = "TCR",
                        align_database = args["align-database"],
                        detailed = args["detailed"],
                        count_chimeric_segments = args["count-chimeric-segments"],
                        chimeric_alignments = args["chimeric-alignments"],
                        )
end

# collect results like so
# igdiscover batch collect --file final/CHMMAIRRa_out_22-1-25.tsv.gz ../  collected_CHMMAIRRa_out.tsv
# pigz -p 20 collected_CHMMAIRRa_out.tsv
