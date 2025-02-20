using CSV, DataFrames, FASTX, Compose, Colors, MolecularEvolution, StringDistances, StringAlgorithms, Distributions, CodecZlib

function read_fasta(filepath::String)
    reader = FASTX.FASTA.Reader(open(filepath, "r"))
    fasta_in = [record for record in reader]
    close(reader)
    return [String(FASTX.FASTA.description(rec)) for rec in fasta_in],
    [uppercase(String(FASTX.FASTA.sequence(rec))) for rec in fasta_in]
end

function write_fasta(filepath::String, sequences::Vector{String}; seq_names = nothing)
    if seq_names === nothing
        seq_names = ["S$(i)" for i in 1:length(sequences)]
    end
    writer = FASTX.FASTA.Writer(open(filepath, "w"))
    for i in 1:length(seq_names)
        rec = FASTX.FASTA.Record(seq_names[i], sequences[i])
        write(writer, rec)
    end
    close(writer)
end

function degap(s::String)
    return replace(s, "-" => "")
end

function parse_recombinations_string(recomb_str::Union{String, Missing})
    if ismissing(recomb_str)
        return []
    end
    recomb_arr = string.(split(recomb_str, ";"))
    return parse_recombination_string.(recomb_arr)
end
function parse_recombination_string(recomb::Union{String, Missing})
    if ismissing(recomb)
        return missing
    end
    recomb = recomb[2:end - 1]
    recomb_split = split(recomb, ", ")
    return (left_allele = recomb_split[1], right_allele = recomb_split[2], left_pos_degapped = parse(Int, recomb_split[3]), right_pos_degapped = parse(Int, recomb_split[4]))
end

function fasttreedbl_nuc_wrapper(seqs::Vector{String}, seqnames::Vector{String}; quiet::Bool = true, binary_name::String = "FastTreeDbl")
    if any(occursin.(":", seqnames))
        @warn "Sanitizing seqnames for FastTreeDbl"
        @warn "Replacing ':' with '_' in seqnames"
        seqnames = [replace(seqname, ":" => "_") for seqname in seqnames]
    end
    mktempdir() do mydir
        seqfile = joinpath(mydir, "sequences.fasta")
        treefile = joinpath(mydir, "fasttree.newick")
        write_fasta(seqfile, seqs, seq_names = seqnames)
        if quiet
            run(`$(binary_name) -quiet -nosupport -nt -gtr -out $(treefile) $(seqfile)`)
        else
            run(`$(binary_name) -nosupport -nt -gtr -out $(treefile) $(seqfile)`)
        end
        lines = open(treefile, "r") do io
            readlines(io)
        end
        return lines[1]
   end
end

function mafft_wrapper(seqs::Vector{String}, seqnames::Vector{String}; mafft::Union{String, Nothing} = nothing, threads::Int = Base.Threads.nthreads())
    # Find mafft executable
    mafft_path = if isnothing(mafft)
        Sys.which("mafft")
    else
        mafft
    end

    if isnothing(mafft_path)
        error("Could not find mafft executable. Please ensure mafft is installed and in your PATH, or provide the path.")
    end

    mktempdir() do mydir
        input_fasta = joinpath(mydir, "sequences.fasta")
        write_fasta(input_fasta, degap.(seqs), seq_names = seqnames)
        cmd = `$(mafft_path) --thread $(threads) $(input_fasta)`
        @info "Running $(cmd)"
        io = IOBuffer()
        run(pipeline(cmd, stdout=io, stderr=devnull))
        stdo = String(take!(io))
        results = split.(split(stdo, ">")[2:end], "\n", limit = 2)
        names = map(x->string(x[1]), results)
        seqs = [uppercase(replace(result[2], "\n" => "")) for result in results ]
        return names, seqs
    end
end

function sanitize_seqname(seqname::String)
    return replace(seqname, ":" => "_")
end

function get_clone_id_sequences(members::DataFrame, clone_id::Int)
    clone_id_members = members[members.clone_id .== clone_id,:];
    seqnames, sequences, chimeric = clone_id_members[!,"sequence_id"], clone_id_members[!,"VDJ_nt"], clone_id_members[!,"chimeric"];
    seqnames = sanitize_seqname.(seqnames);
    sort!(clone_id_members, :V_SHM)
    root_seq = clone_id_members[.! clone_id_members.chimeric,"germline_alignment"][1]
    sequences = vcat(root_seq, sequences)
    seqnames = vcat("root", seqnames)
    chimeric = vcat(false, chimeric)
    return seqnames, sequences, chimeric
end

# align seqs, get tree, read into tree object, reroot and ladderize
function seqs2fasttree_tree(sequences::Vector{String}, seqnames::Vector{String})
    _, aligned_sequences = mafft_wrapper(sequences, seqnames)
    str_tree = fasttreedbl_nuc_wrapper(aligned_sequences, seqnames)
    newt = gettreefromnewick(str_tree, FelNode)
    newt = MolecularEvolution.reroot!([x for x in getnodelist(newt) if x.name=="root"][1], dist_above_child=0)
    ladderize!(newt)
    return newt
end

function countmap_df(arr)
    d = DataFrame([(el[1], el[2]) for el in countmap(arr)])
    rename!(d, :1 => "value", :2 => "n")
    sort!(d, :n, rev = true)
    return d
end

function fasta2df(fastapath)
    n,s = read_fasta(fastapath)
    return DataFrame(:sequence_id => n, :sequence => s)
end

function inverse_longestcommonsubstring(s1::String, s2::String)
    return 1 / length(longestcommonsubstring(s1, s2)[1])
end

function lev(s1, s2)
    return Levenshtein()(s1, s2)
end

# takes in alleles and a distance metric
# produces a gene-averaged hclust result
function alleles2distmats(sequence_ids::Vector, sequences::Vector; dist = inverse_longestcommonsubstring)
    distmat = pairwise(dist, sequences, sequences)
    # collapse the pairwise allele distmat into a pairwise gene distmat
    genes = map(x->split(x, "*")[1], sequence_ids)
    unique_genes = sort(unique(genes))
    gene_mat = zeros(length(unique_genes), length(unique_genes))
    for i in 1:length(unique_genes)
        for j in i+1:length(unique_genes)
            inds1 = findall(x->x == unique_genes[i], genes)
            inds2 = findall(x->x == unique_genes[j], genes)
            coords = [[ind1, ind2] for ind1 in inds1 for ind2 in inds2]
            gene_mat[j, i] = gene_mat[i, j] = mean([distmat[coord[1], coord[2]] for coord in coords])
        end
    end
    return (allele_distances = distmat, gene_distances = gene_mat, genes = unique_genes, alleles = sequence_ids)
end

function alleles2genehclust(sequence_ids::Vector, sequences::Vector; dist = inverse_longestcommonsubstring)
    distmats = alleles2distmats(sequence_ids, sequences, dist = dist)
    h = hclust(distmats.gene_distances, linkage = :average, branchorder = :optimal)
    return (hclust = h, hclust_genes = distmats.genes[h.order], gene_distances = distmats.gene_distances[h.order,h.order], allele_distances = distmats.allele_distances[h.order,h.order])
end

function upperright_sum(mat)
    s = 0
    for i in 1:size(mat)[1]
        for j in i+1:size(mat)[1]
            s += mat[i,j]
        end
    end
    return s
end

function lowerleft_sum(mat)
    s = 0
    for i in 1:size(mat)[1]
        for j in 1:i - 1
            s += mat[i,j]
        end
    end
    return s
end

function upperright_vec(mat)
    s = 0
    v = Float64[]
    for i in 1:size(mat)[1]
        v = vcat(v, mat[i,i+1:size(mat)[1]])
    end
    return v
end

function lowerleft_vec(mat)
    s = 0
    v = Float64[]
    for i in 1:size(mat)[1]
        v = vcat(v, mat[i,1:i-1])
    end
    return v
end

function add_jitter(num; σ = 0.05, μ = 0)
    return num + rand(Normal(μ, σ))
end

function jitter(a::Array; variance=1.0) 
    return a .+ rand(truncated(Normal(0, variance), -.1, 0.1), size(a)) 
end



function inverse_longestcommonsubstring(s1::String, s2::String)
    return 1 / length(longestcommonsubstring(s1, s2)[1])
end

# library-related augmentation
function call2gene(call) return ismissing(call) ? missing : split(call, "*")[1] end
function gene2family(gene) return ismissing(gene) ? missing : split(gene, "-")[1] end
function library2chain(library)
    pattern2chain = Dict("TRA" => "TRA", "TRB" => "TRB", "TRD" => "TRD", "TRG" => "TRG", "IgG" => "IGG", "IgM" => "IGM", "IGM" => "IGM", "IGG" => "IGG")
    for pair in collect(pattern2chain)
        if occursin(pair[1], library)
            return pair[2]
        end
    end
    return ""
end

"""
Read in dataset collected files, concatenating and adding dataset_id
"""
function load_collected_files(collected_files)
    chimeric_concat_df = DataFrame()
    nonchimeric_concat_df = DataFrame()
    library_cts_concat_df = DataFrame()
    for collected_file in collected_files
        df = DataFrame(CSV.File(collected_file, delim = '\t'; drop = ["threaded"]))
        dataset_id = split(collected_file, "/")[end-2]
        @info dataset_id
        df[!,"dataset_id"] .= dataset_id
        df[!,"chimeric"] .= df.chimera_probability .> 0.95
        append!(chimeric_concat_df, df[df.chimeric, :], cols = :union)
        # we only want the startingpoint values from the nonchimeras
        # loading the whole thing takes too much memory
        nonchimeric_df = combine(groupby(df[.!df.chimeric, :], [:dataset_id, :case, :startingpoint]), :case => length => :n)
        append!(nonchimeric_concat_df, nonchimeric_df, cols = :union)
        library_cts = combine(groupby(df, :case), :case => length => :library_reads, :dataset_id => first => :dataset_id)
        append!(library_cts_concat_df, library_cts, cols = :union)
    end
    return chimeric_concat_df, nonchimeric_concat_df, library_cts_concat_df
end


"""
Add columns to collected files
"""
function augment_collected_df(collected)
    collected[!,"locus"] .= map(x->x[1:3], collected[!,"startingpoint"]);
    collected[!,"v_gene"] .= call2gene.(collected[!,"startingpoint"]);
    collected[!,"v_family"] .= gene2family.(collected[!,"v_gene"]);
    collected[!,"chain"] = library2chain.(collected[!,"case"]);
    if "dataset_id" in names(collected)
        collected[!,"donor"] .= "unknown"
        collected[collected.dataset_id .== "GKH_TCR","donor"] .= map(x->split(x, "_")[2], collected[collected.dataset_id .== "GKH_TCR","case"])
    end
    function jitter(a::Array; variance=1.0) return a .+ rand(truncated(Normal(0, variance), -.1, 0.1), size(a)) end
    collected[!,"chain_coord"] = jitter(map(function(x) chain2coord[x] end, collected[!,"chain"]), variance = 2);
    return collected
end


function get_igdiscover_run_dirs(d::String)::Array{String,1} return [abspath(el[1]) for el in walkdir(d, follow_symlinks = true) if "igdiscover.yaml" in el[3]] end
function numlines(f::String)::Int64
    if isfile(f)
        return f[end - 2:end] == ".gz" ?  count(x->x == '\n', String(read(`zcat $f`))) : countlines(f)
    end
    return 0
end

function gather_igdiscover_counts(igdiscover_dataset_dir::String, chimeric_file::String; p_threshold::Float64 = 0.95)::DataFrame
    """
    Takes a directory containing a set of IgDiscover runs
    Returns a dataframe with the library name, total reads, merged reads, assigned reads, filtered reads.
    """
    println(igdiscover_dataset_dir)
    println(chimeric_file)
    dataset = basename(igdiscover_dataset_dir)
    read_cts_df = DataFrame(dataset = String[], library = String[], total = Int[], merged = Int[], assigned = Int[], filtered = Int[], chimeric = Int[], percent_chimeric = Float16[])
    igdiscover_dirs = get_igdiscover_run_dirs(igdiscover_dataset_dir)
    cases = [join(split(relpath(igdiscover_dir, igdiscover_dataset_dir), "/"), "_") for igdiscover_dir in igdiscover_dirs]
    for (case, igdiscover_dir) in zip(cases, igdiscover_dirs)
        f = joinpath(igdiscover_dir, "stats/reads.json")
        if ! isfile(f)
            @warn "No such file $(f)"
            continue
        end
        reads_JSON = JSON.parsefile(f)
        delete!(reads_JSON, "merging_was_done")
        reads_JSON["library"] = case
        assigned_JSON = JSON.parsefile(joinpath(igdiscover_dir, "final/stats/assigned.json"))
        reads_JSON["assigned"] = assigned_JSON["total"]
        filtered_JSON = JSON.parsefile(joinpath(igdiscover_dir, "final/stats/filtered.json"))
        reads_JSON["filtered"] = filtered_JSON["has_cdr3"]
        reads_JSON["dataset"] = dataset
        chimera_probabilities = CSV.read(joinpath(igdiscover_dir, chimeric_file), DataFrame, delim = '\t', select = ["chimera_probability"])[!,"chimera_probability"]
        reads_JSON["chimeric"] = sum(chimera_probabilities .>= p_threshold)
        reads_JSON["percent_chimeric"] = reads_JSON["chimeric"]  / reads_JSON["filtered"] .* 100
        push!(read_cts_df, reads_JSON)
    end
    read_cts_df[!,:chain] = library2chain.(read_cts_df.library)
    return read_cts_df
end

function library2chain(library::String)::String
    chains = ["TRA", "TRB", "TRD", "TRG", "IgG", "IgM", "IGM"]
    for chain in chains
        if occursin(chain, library)
            return chain
        end
    end
    return ""
end

function run_igblastwrap_from_files(db_dir, fastq_path, out_path; sequence_type = "IG", human_gl_aux = nothing, threads = Base.Threads.nthreads())
    output = IOBuffer()
    run(pipeline(`conda run -n igdiscover igdiscover igblastwrap --sequence-type $(sequence_type) --aux $(human_gl_aux) --threads $(threads) $(db_dir) $(fastq_path)`, output))
    if out_path[end - 1 : end] == "gz"
        open(GzipCompressorStream, out_path, "w") do io
            write(io, take!(output))
        end
    else
        open(out_path, "w") do io
            write(io, take!(output))
        end
    end
end

function run_augment_from_files(db_dir, assignments_path, out_path)
    output = IOBuffer()
    run(pipeline(`conda run -n igdiscover igdiscover augment --read-cdr3 $(db_dir) $(assignments_path)`, output))
    if out_path[end - 1 : end] == "gz"
        open(GzipCompressorStream, out_path, "w") do io
            write(io, take!(output))
        end
    else
        open(out_path, "w") do io
            write(io, take!(output))
        end
    end
end

function run_filter_from_files(assignments_path, out_path)
    output = IOBuffer()
    run(pipeline(`conda run -n igdiscover igdiscover filter $(assignments_path)`, output))
    if out_path[end - 1 : end] == "gz"
        open(GzipCompressorStream, out_path, "w") do io
            write(io, take!(output))
        end
    else
        open(out_path, "w") do io
            write(io, take!(output))
        end
    end
end