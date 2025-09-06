using CSV, DataFrames, Random, StatsBase, Plots

# samples n alleles per gene from the reference set
function simulate_genotype(refnames, refseqs, rng; n_alleles_per_gene = 1)
    genenames = map(refname -> split(refname, "*")[1], refnames)
    unique_genenames = unique(genenames)
    u_gene_inds = vcat([sample(rng, findall(genenames .== genename), n_alleles_per_gene) for genename in unique_genenames]...)
    # take one allele per gene
    refnames, refseqs = refnames[u_gene_inds], refseqs[u_gene_inds]
    
    refnames, refseqs = mafft_wrapper(refseqs, refnames, parameters = ["--localpair"])
    return refnames, refseqs
end

function random_seqpair(names::Vector{String}, seqs::Vector{String}, rng)
    rand_inds = sample(rng, 1:length(names), 2, replace = false)
    return names[rand_inds], seqs[rand_inds]
end

function get_breakpoint(seqlen::Int, rng; min_pos::Union{Int, Float64} = 0.0, max_pos::Union{Int, Float64} = 1.0)
    if ismissing(max_pos)
        max_pos = seqlen
    end
    if (typeof(min_pos) == Float64) & (min_pos < 1)
        min_pos = floor(min_pos * seqlen)
    end
    if (typeof(max_pos) == Float64) & (max_pos <= 1)
        max_pos = floor(max_pos * seqlen)
    end
    if min_pos == max_pos
        return Int(min_pos)
    end
    return Int(rand(rng, min_pos:max_pos))
end

function add_random_shm(seq::String, rng; min_shm::Float64 = 0.01, max_shm::Float64 = 0.2)
    if (min_shm > 1.0) | (max_shm > 1.0) | (min_shm < 0.0) | (max_shm < 0.0) | (min_shm > max_shm)
        error("min_shm and max_shm must be between 0 and 1, and min_shm must be less than max_shm")
    end
    NUCLEOTIDES = ['A', 'C', 'G', 'T']
    seq = collect(join(seq))
    # get a count of shm positions based on the shm rate
    shm_rate = rand(rng, min_shm:0.0001:max_shm)
    shm_ct = Int(floor(length(seq) * shm_rate))
    # choose positions to apply shm to, exclude gaps
    shm_positions = shuffle([i for i in 1:length(seq) if seq[i] != '-'])[1:shm_ct]
    # apply shm by changing to a random nucleotide that is not the current nucleotide
    for shm_position in shm_positions
        curr_nuc = seq[shm_position]
        seq[shm_position] = setdiff(NUCLEOTIDES, curr_nuc)[randperm(3)[1]]
    end
    return join(seq)
end


"""
Adds illumina MiSeq errors to a set of sequences in 3 steps: 
1. Add prefix and suffix to each sequence to make it the length of the real template amplified by the Illumina sequencer
2. Generate paired end reads with errors using art_illumina
3. Merge the two reads into a single read using pear

The output is a vector of tuples, each containing the sequence ID, sequence, and quality score for each read.

Note:
- specifically for Illumina MiSeq v3
- generates 250bp reads
"""
function add_illumina_errors(seqs::Vector{String}, prefix_nucs::Int64, suffix_nucs::Int64; art_illumina::Union{String, Nothing} = nothing, pear_path::Union{String, Nothing} = nothing, gene::Char = 'v', padding_dists = Dict("before_V" => 40, 'v' => 290, 'd' => 20, 'j' => 50, "after_J" => 40), V_seq = "", D_seq = "", J_seq = "")

    art_illumina_path = if isnothing(art_illumina)
        Sys.which("art_illumina")
    else
        art_illumina
    end

    if isnothing(art_illumina_path)
        error("Could not find art_illumina executable. Please ensure art_illumina is installed and in your PATH, or provide the path.")
    end

    seqs = degap.(seqs)
    prefix_string, suffix_string = "", ""
    # 1. Add prefix and suffix
    # we need to pad the sequences to place them in the correct positions in the amplicon for realistic errßor simulation
    if gene == 'v'
        prefix_string = join([rand_nuc() for i in 1:padding_dists["before_V"]])
        suffix_string = string(D_seq, J_seq, join([rand_nuc() for i in 1:padding_dists["after_J"]]))
    elseif gene == 'd'
        prefix_string = string(join([rand_nuc() for i in 1:padding_dists["before_V"]]), V_seq)
        suffix_string = string(J_seq, join([rand_nuc() for i in 1:padding_dists["after_J"]]))
    elseif gene == 'j'
        prefix_string = string(join([rand_nuc() for i in 1:padding_dists["before_V"]]), V_seq, D_seq)
        suffix_string = join([rand_nuc() for i in 1:padding_dists["after_J"]])
    else
        error("Invalid gene: $(gene)")
    end

    seqs = [string(prefix_string, seq, suffix_string) for seq in seqs]
    total_prefix_length = length(prefix_string)
    total_suffix_length = length(suffix_string)
    

    mktempdir() do mydir
        # 2. Generate paired end reads with errors
        input_fasta = joinpath(mydir, "sequences.fasta")
        write_fasta(input_fasta, seqs, seq_names = ["$(i)" for i in 1:length(seqs)])
        cmd = `$(art_illumina_path) -k 0 -ss MSv3 -amp -p -i $(input_fasta) -o $(joinpath(mydir, "amplicon")) -f 1 -l 250`
        io = IOBuffer()
        run(pipeline(cmd, stdout=io, stderr=devnull))
        # 3. Merge generated reads
        merged_sequence_ids, merged_sequences, merged_qualities = read_fastq(merge_fastqs_from_files("$(mydir)/amplicon1.fq", "$(mydir)/amplicon2.fq", "$(mydir)/merged", pear = pear_path))
        merged_sequences = [seq[total_prefix_length + 1 : end - total_suffix_length] for seq in merged_sequences]
        return merged_sequences
    end
end

function merge_fastqs_from_files(r1_filepath::String, r2_filepath::String, merged_filepath::String; pear::Union{String, Nothing} = nothing)
    pear_path = if isnothing(pear)
        Sys.which("pear")
    else
        pear
    end

    if isnothing(pear_path)
        error("Could not find pear executable. Please ensure pear is installed and in your PATH, or provide the path.")
    end
    run(Cmd(String.([split(pear_path, " ")..., "-f", r1_filepath, "-r", r2_filepath, "-o", merged_filepath])))
    return merged_filepath * ".assembled.fastq"
end


function add_variations(seqs::Vector{String}, rng; min_shm::Float64 = 0.01, max_shm::Float64 = 0.2, variation_method::String = "uniform_random", shazam_implementation::String = "parallel", prefix_nucs::Int64 = 40, suffix_nucs::Int64 = 80, gene::Char = 'v', V_seq = "", D_seq = "", J_seq = "")
    if variation_method == "uniform_random"
        seqs = map(seq -> add_random_shm(seq, rng, min_shm = min_shm, max_shm = max_shm), seqs) # add shm
    elseif variation_method == "shazam"
        if min_shm != max_shm
            error("min_shm and max_shm must be the same for shazam shmulateSeq")
        end
        # Choose shazam implementation
        if shazam_implementation == "original"
            seqs = add_shazam_shm(seqs, min_shm)
        elseif shazam_implementation == "vectorized"
            seqs = add_shazam_shm_vectorized(seqs, min_shm)
        elseif shazam_implementation == "parallel"
            seqs = add_shazam_shm_parallel(seqs, min_shm)
        else
            error("Invalid shazam_implementation: $(shazam_implementation). Use 'original', 'vectorized', or 'parallel'")
        end
    elseif variation_method == "art_illumina"
        seqs = add_illumina_errors(seqs, prefix_nucs, suffix_nucs, pear_path = "conda run -n igdiscover pear", gene = gene, V_seq = V_seq, D_seq = D_seq, J_seq = J_seq)
    else
        error("Invalid variation_method: $(variation_method)")
    end
    return seqs
end


function add_shazam_shm_parallel(seqs::Vector{String}, shm_rate::Float64; cores::Int = 5)
    # if the sequence is aligned, remove those gaps, apply shazam, then add the gaps back
    nongap_positions = findall.(x->x != '-', seqs)
    degapped_seqs = degap.(seqs)

    res = mktempdir() do dir
        open("$(dir)/seqs.txt", "w") do f
            write(f, join(degapped_seqs, "\n"))
        end
        r_string = """
        library(shazam)
        library(parallel)
        arr = readLines("$(dir)/seqs.txt")
        # Use parallel processing for large datasets
        cl <- makeCluster($(cores), type = "FORK")
        clusterEvalQ(cl, library(shazam))
        new_arr = parSapply(cl, arr, function(seq) shmulateSeq(seq, $(shm_rate), frequency = TRUE))
        stopCluster(cl)
        writeLines(paste(new_arr, collapse = "\n"), "$(dir)/output.txt")"""
        run(`Rscript -e $r_string`) 
        return read("$(dir)/output.txt", String)
    end
    res = String.(split(strip(res), "\n"))

    for i in 1:length(seqs)
        seq = collect(seqs[i])
        seq[nongap_positions[i]] = collect(res[i])
        seqs[i] = String(seq)
    end

    return seqs
end

function add_shazam_shm_vectorized(seqs::Vector{String}, shm_rate::Float64)
    seqs = [replace(seq, '-' => '.') for seq in seqs]

    local res = ""
    mktempdir() do dir
        open("$(dir)/seqs.txt", "w") do f
            write(f, join(seqs, "\n"))
        end
        r_string = """
        library(shazam)
        arr = readLines("$(dir)/seqs.txt")
        # Use sapply for vectorization - much faster than explicit loops
        new_arr = sapply(arr, function(seq) shmulateSeq(seq, $(shm_rate), frequency = TRUE), USE.NAMES = FALSE)
        writeLines(paste(new_arr, collapse = "\n"), "$(dir)/output.txt")"""
        run(`Rscript -e $r_string`) 
        res = read("$(dir)/output.txt", String)
    end
    res = [replace(seq, '.' => '-') for seq in String.(split(strip(res), "\n"))]
    return res
end

function add_shazam_shm(seqs::Vector{String}, shm_rate::Float64)
    seqs = [replace(seq, '-' => '.') for seq in seqs]

    local res = ""
    mktempdir() do dir
        open("$(dir)/seqs.txt", "w") do f
            write(f, join(seqs, "\n"))
        end
        r_string = """
        library(shazam)
        arr = readLines("$(dir)/seqs.txt")
        new_arr = character(length(arr))  # Pre-allocate vector
        for (i in 1:length(arr)){
            new_arr[i] = shmulateSeq(arr[i], $(shm_rate), frequency = TRUE)
        }
        writeLines(paste(new_arr, collapse = "\n"), "$(dir)/output.txt")"""
        run(`Rscript -e $r_string`) 
        res = read("$(dir)/output.txt", String)
    end
    res = [replace(seq, '.' => '-') for seq in String.(split(strip(res), "\n"))]
    return res
end

# pick n random sequences from the reference set and add shm to them
# shm rate chosen uniformly between min_shm and max_shm
# assumes sequences are already aligned
function random_nonchimeras(names::Vector{String}, seqs::Vector{String}, rng; n::Int = 1, min_shm::Float64 = 0.2, max_shm::Float64 = 0.2, variation_method::String = "uniform_random", gene::Char = 'v', V_seq = "", D_seq = "", J_seq = "", pre_aligned::Bool = true)
    #if length(unique(length.(seqs))) != 1
    #    error("All sequences must be the same length. Please align them and try again.")
    #end
    rand_inds = sample(rng, 1:length(names), n, replace = true)
    rand_names, rand_seqs = names[rand_inds], seqs[rand_inds]
    rand_seqs = add_variations(rand_seqs, rng, min_shm = min_shm, max_shm = max_shm, variation_method = variation_method, gene = gene, V_seq = V_seq, D_seq = D_seq, J_seq = J_seq)
    rand_names = ["$(rand_names[i])_$(i)" for i in 1:length(rand_names)]
    return rand_names, rand_seqs
end

# Generate n chimeras by taking random sequences, cutting at a random breakpoint between min_pos and max_pos, and stitching together the cut sequences
# assumes sequences are already aligned
function random_chimeras(names::Vector{String}, seqs::Vector{String}, rng; min_pos::Union{Int, Float64} = 0.0, max_pos::Union{Int, Float64} = 1.0, n::Int = 1, min_shm1::Float64 = 0.0, max_shm1::Float64 = 0.2, min_shm2::Float64 = 0.0, max_shm2::Float64 = 0.2, variation_method::String = "uniform_random", V_seq = "", D_seq = "", J_seq = "", gene::Char = 'v', ig_seqtype = "TCR", pre_aligned::Bool = true)
    # we generate 2x as many sequences so we can throw away the ones that don't align well with IgBLAST and then take n
    preliminary_n = n * 100
    rand_refs_list = [random_seqpair(names, seqs, rng) for i in 1:preliminary_n] # choose seqeunces
    seqs1 = add_variations(map(rand_refs -> rand_refs[2][1], rand_refs_list), rng, min_shm = min_shm1, max_shm = max_shm1, variation_method = variation_method, gene = gene, V_seq = V_seq, D_seq = D_seq, J_seq = J_seq)
    seqs2 = add_variations(map(rand_refs -> rand_refs[2][2], rand_refs_list), rng, min_shm = min_shm2, max_shm = max_shm2, variation_method = variation_method, gene = gene, V_seq = V_seq, D_seq = D_seq, J_seq = J_seq)
    rand_refs_list = [[rand_refs_list[i][1], [seqs1[i], seqs2[i]]] for i in 1:length(rand_refs_list)] 

    recombined = [recombine_seqs(rand_refs[1][1], rand_refs[1][2], rand_refs[2][1], rand_refs[2][2], rng, min_pos, max_pos, pre_aligned = pre_aligned) for rand_refs in rand_refs_list]
    recombined_names, recombined_seqs, recombined_breakpoints = map(x-> x[2], recombined), map(x-> x[1], recombined), map(x-> x[3], recombined)
    recombined_names = ["$(recombined_names[i])_$(i)" for i in 1:length(recombined_names)]

    if ! pre_aligned
        # run igblast so we can exclude chimeras that don't align to the reference set
        assignments = DataFrame()
        if gene == 'v'
            assignments = igblast_v(recombined_seqs, recombined_names, degap.(seqs), namesCS, D_seq = D_seq, J_seq = J_seq, ig_seqtype = ig_seqtype)
        elseif gene == 'd'
            assignments = igblast_d(recombined_seqs, recombined_names, degap.(seqs), names, V_seq = V_seq, J_seq = J_seq, ig_seqtype = ig_seqtype)
        elseif gene == 'j'
            assignments = igblast_j(recombined_seqs, recombined_names, degap.(seqs), names, V_seq = V_seq, D_seq = D_seq, ig_seqtype = ig_seqtype)
        end
        assignments = assignments[ .! ismissing.(assignments[!, "$(gene)_call"]),:]
        assignments[!, "$(gene)_covered"] .= [covered(record["$(gene)_germline_alignment"], seqs[names .== record["$(gene)_call"]][1]) for record in eachrow(assignments)]
        covered_mask = assignments[!, "$(gene)_covered"] .> 90.0
        chimeric_names, chimeric_seqs, breakpoint_positions = recombined_names[covered_mask][1:n], recombined_seqs[covered_mask][1:n], recombined_breakpoints[covered_mask][1:n]
        return chimeric_names, chimeric_seqs, breakpoint_positions
    else
        return recombined_names, recombined_seqs, recombined_breakpoints
    end
end

function recombine_seqs(name1::String, name2::String, seq1::String, seq2::String, rng, min_pos, max_pos; pre_aligned::Bool = true)
    if ! pre_aligned
        scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1)
        res = pairalign(GlobalAlignment(), seq1, seq2, scoremodel)
        aligned_seq1 = res.aln.a.seq
        aligned_seq2 = res.aln.b
        breakpoint = get_breakpoint(length(aligned_seq1), rng, min_pos = min_pos, max_pos = max_pos)
        recombined_seq = degap(string(aligned_seq1[1:breakpoint], aligned_seq2[breakpoint + 1:end]))
        recombined_name = "$(name1)_1to$(breakpoint)_$(name2)_$(breakpoint)toEnd" 
    else
        breakpoint = get_breakpoint(length(seq1), rng, min_pos = min_pos, max_pos = max_pos)
        recombined_seq = string(seq1[1:breakpoint], seq2[breakpoint + 1:end])
        recombined_name = "$(name1)_1to$(breakpoint)_$(name2)_$(breakpoint)toEnd" 
    end
    return recombined_seq, recombined_name, breakpoint
end



function usearch_uchime2_ref_wrapper(query_names::Vector{String}, query_seqs::Vector{String}, db_seqs::Vector{String}; search::String = "usearch", mode = "sensitive", mindiv = 0.00001)
    """
    Runs uchime_ref and returns sequence scores
    Note 1: query_names must be unique because we need them to sort the uchime_ref output
    Note 2: I was running usearch using a docker container because usearch doesn't work on apple silicon. You may be able to rewrite the command in a simpler way. Binaries are over here https://drive5.com/usearch/download.html
    """
    local res
    mktempdir() do dir
        # run uchime_ref
        db_fasta_path = joinpath(dir, "V.fasta")
        query_fasta_path = joinpath(dir, "queries.fasta")
        outpath = joinpath(dir, "out.txt")
        write_fasta(db_fasta_path, degap.(db_seqs))
        write_fasta(query_fasta_path, degap.(query_seqs), seq_names = query_names)
        #cmd = `usearch -uchime2_ref $(query_fasta_path) -db $(db_fasta_path) -uchimeout $(outpath) -strand plus -mode $(mode)`
        cmd = Cmd([
            "docker", "run", "--rm", "--interactive", "--platform", "linux/amd64",
            "--volume", "$(dir):/workdir", "--env", "ubuntu", "multiarch/crossbuild",
            "/bin/bash", "-c", """
            cd /workdir
            # grab usearch binary
            wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz
            gunzip usearch11.0.667_i86linux32.gz
            chmod +x usearch11.0.667_i86linux32

            # run usearch
            ./usearch11.0.667_i86linux32 -uchime2_ref queries.fasta -db V.fasta -uchimeout out.txt -strand plus -mode $(mode) -mindiv $(mindiv)
            """
        ])
        
        
        println(cmd)
        run(cmd)
        # get results and align them with the labels
        res = CSV.read(outpath, delim = "\t", DataFrame, header = [ "sequence_id", "score", "chimera", "chimera_from", "chimera_to", "model", "attributes"])
    end
    order_lookup = Dict(zip(query_names, collect(1:length(res.sequence_id))))
    res = res[sortperm(res.sequence_id, by=x->order_lookup[x]),:]
    attribute_columns = vcat(parse_attributes.(res.attributes)...)
    res = hcat(res, attribute_columns)
    return res
end

function parse_attributes(attributes::AbstractString)
    attribute_dict::Dict{String, Any} = Dict("dqt" => missing, "dqm" => missing, "div" => missing, "L_Y" => missing, "L_N" => missing, "L_A" => missing, "R_Y" => missing, "R_N" => missing, "R_A" => missing, "L" => missing, "R" => missing, "why" => missing)
    if occursin(";", attributes)
        for attribute in [el for el in split(attributes, ";") if el != ""]
            n,v = split(attribute, "=")
            if n == "dqt"
                v = parse(Int, v)
            elseif n == "dqm"
                v = parse(Int, v)
            elseif n == "div"
                v = parse(Float64, v[1:end - 1])
            elseif (n == "L") | (n == "R")
                attribute_dict["$(n)_Y"] = parse(Int, split(v, (',', '('))[1])
                attribute_dict["$(n)_N"] = parse(Int, split(v, (',', '('))[2])
                attribute_dict["$(n)_A"] = parse(Int, split(v, (',', '('))[3])
            end
            attribute_dict[n] = v
        end
    end
    return DataFrame(attribute_dict)
end

function vsearch_uchime_ref_wrapper(query_names::Vector{String}, query_seqs::Vector{String}, db_seqs::Vector{String}; vsearch::String = "vsearch", mindiv::Float64 = 0.8, mindiffs::Int = 3, xn::Int = 8)
    """
    Runs uchime_ref and returns sequence scores
    Note that query_names must be unique because we need them to sort the uchime_ref output
    """
    @info "Running vsearch_uchime_ref_wrapper"
    local res
    mktempdir() do dir
        # run uchime_ref
        uchime_nonchimeras_path = joinpath(dir, "uchime_ref_nonchimeras.fasta")
        uchime_chimeras_path = joinpath(dir, "uchime_ref_chimeras.fasta")
        db_fasta_path = joinpath(dir, "V.fasta")
        query_fasta_path = joinpath(dir, "queries.fasta")
        outpath = joinpath(dir, "out.txt")
        write_fasta(db_fasta_path, degap.(db_seqs))
        write_fasta(query_fasta_path, degap.(query_seqs), seq_names = query_names)
        cmd = `vsearch --uchime_ref $(query_fasta_path) --uchimeout $(outpath) --fasta_score --db $(db_fasta_path) --mindiv $(mindiv) --mindiffs $(mindiffs) --xn $(xn)`
        println(cmd)
        run(cmd)
        cp(outpath, "/Users/march712/Downloads/vsearch_uchime_out.txt", force = true)
        # get results and align them with the labels
        res = CSV.read(outpath, delim = "\t", DataFrame, header = ["score", "sequence_id", "parent_A", "parent_B", "top_parent", "idQM", "idQA", "idQB", "idAB", "idQT", "LY", "LN", "LA", "RY", "RN", "RA", "div", "YN"])
    end
    order_lookup = Dict(zip(query_names, collect(1:length(res.sequence_id))))
    scores = res.score
    scores = scores[sortperm(res.sequence_id, by=x->order_lookup[x])]
    return scores
end

function get_confusion_tuple(results, labels)
    return (TP = sum(results .& labels), FP = sum(results .& (.! labels)), FN = sum((.! results) .& labels), TN = sum((.! results) .& (.! labels)))
end

function get_rates(results, labels)
    conf = get_confusion_tuple(results, labels)
    return (TPR = conf.TP / (conf.TP + conf.FN), FNR = conf.FN / (conf.FN + conf.TP), TNR = conf.TN / (conf.TN + conf.FP), FPR = conf.FP / (conf.FP + conf.TN))
end

function usearch_uchime2_ref_ROC_curve(query_names::Vector{String}, query_seqs::Vector{String}, labels, ref_seqs::Vector{String}; mindiffs = 2, mindivt = 0.5, xa = 1, xn = 4, modes = ["high_confidence", "specific", "balanced", "sensitive"], mindiv = 0.00001)
    """
    Generate an ROC curve plot for a set of query sequences using the CHMMera method by varying the chimeric probability threshold
    """
    calls = Dict()
    # we'll want the last output from the loop, so we need to initialize the variables outside of the loop
    TPRs, FPRs, scores, cutoffs = zeros(0), zeros(0), zeros(0), zeros(0)
    for mode in modes
        res = usearch_uchime2_ref_wrapper(query_names, query_seqs, ref_seqs, mode = mode, mindiv = usearch_uchime2_mindiv)
        scores, chimera = res.score, res.chimera
        # gather TPRs and FPRs for different cutoffs
        cutoffs = collect(0:.01:maximum(scores))
        TPRs, FPRs = zeros(0), zeros(0)
        for prob_cutoff in cutoffs
            no_missing = (.! ismissing.(res.div)) .& (.! ismissing.(res.L_Y)) .& (.! ismissing.(res.R_Y)) .& (.! ismissing.(res.score))
            results = no_missing .& (res.score .> prob_cutoff) .& (res.L_Y .>= mindiffs) .& (res.R_Y .>= mindiffs) .& (res.div .>= mindivt)
            rates = get_rates(results, labels)
            push!(TPRs, rates.TPR)
            push!(FPRs, rates.FPR)
        end
        rates = get_rates(chimera .== "Y", labels)
        calls[mode] = (TPR = rates.TPR, FPR = rates.FPR)
    end
    println("MINIMUM USEARCH SCORE $(minimum(scores))")
    # the built-in cutoffs in USEARCH prevent us from reaching TPR/FPR 1, so I add it manually
    prepend!(FPRs, 1.0)
    prepend!(TPRs, 1.0)
    push!(scores, 0.000000001)
    push!(cutoffs, 0.000000001)
    return FPRs, TPRs, scores, cutoffs, calls
end

function vsearch_uchime_ref_ROC_curve(query_names::Vector{String}, query_seqs::Vector{String}, labels, ref_seqs::Vector{String}; cutoff_interval = 0.01, mindiv = 0.8, mindiffs = 3, xn = 8)
    """
    Generate an ROC curve plot for a set of query sequences using the CHMMera method by varying the chimeric probability threshold
    """
    # get uchime_ref scores
    scores = vsearch_uchime_ref_wrapper(query_names, query_seqs, ref_seqs, mindiv = vsearch_uchime_mindiv, mindiffs = vsearch_uchime_mindiffs)
    # gather TPRs and FPRs for different cutoffs
    @info "DONE WITH VSEARCH UCHIME"
    cutoffs = collect(0:cutoff_interval:maximum(scores))
    rates = map(cutoff -> get_rates(scores .>= cutoff, labels), cutoffs)
    @info "DONE WITH RATES"
    return map(x-> x.FPR, rates), map(x-> x.TPR, rates), scores, cutoffs
end

function CHMMera_ROC_curve(query_seqs::Vector{String}, labels, ref_seqs::Vector{String}; bw = true, mutation_probabilities = [0.001, 0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2], base_mutation_probability = 0.05, prior_probability = 0.02, cutoff_interval = 0.01)
    """
    Generate an ROC curve plot for a set of query sequences using the CHMMera method by varying the chimeric probability threshold
    """
    rates = []
    # new format
    chimera_probs = CHMMera.get_chimera_probabilities(query_seqs, ref_seqs, bw = bw, mutation_probabilities = mutation_probabilities, base_mutation_probability = base_mutation_probability, prior_probability = prior_probability)
    # gather TPRs and FPRs for different cutoffs
    cutoffs = collect(0:cutoff_interval:1.1)
    rates = map(cutoff -> get_rates(chimera_probs .>= cutoff, labels), cutoffs)
    return map(x-> x.FPR, rates), map(x-> x.TPR, rates), chimera_probs, cutoffs
end

function CHMMAIRRa_ROC_curve(assignments::DataFrame, labels, ref_seqs::Vector{String}, ref_names::Vector{String}; cutoff_interval = 0.01, receptor = "TCR", gene::Char = 'v', V_seq = "", D_seq = "", J_seq = "", bw = true, ig_seqtype = "TCR", prior_probability = 0.05)
    mutation_probabilities = Float64[]
    if (receptor == "IG") & (! bw)
        mutation_probabilities = [0.001, 0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2]
    elseif (receptor == "TCR") & ( ! bw)
        mutation_probabilities = [0.005]
    end
    HMM_parameters = Dict("method" => bw ? "BW" : "DB", "mutation_probabilities" => mutation_probabilities, "prior_probability" => prior_probability, "base_mutation_probability" => 0.05)

    res = CHMMAIRRa.detect_chimeras(ref_names, ref_seqs, assignments, gene, receptor = receptor, HMM_parameters = HMM_parameters)

    assignments = leftjoin(assignments, res.out[!,["sequence_id", "$(gene)_chimera_probability", "$(gene)_chimeric", "$(gene)_threaded"]], on = :sequence_id)
    cutoffs = collect(0:cutoff_interval:1.1)
    rates = map(cutoff -> get_rates(assignments[!, "$(gene)_chimera_probability"] .>= cutoff, labels), cutoffs)
    return map(x-> x.FPR, rates), map(x-> x.TPR, rates), assignments[!, "$(gene)_chimera_probability"], cutoffs, assignments
end

function AUC(TPRs::Vector{Float64}, FPRs::Vector{Float64})
    """Calculate AUC using trapezoids"""
    AUC = 0
    for i in 2:length(TPRs)
        AUC = AUC +  ((FPRs[i] - FPRs[i - 1]) * ((TPRs[i] + TPRs[i - 1]) / 2))
    end
    return abs(AUC)
end




# this one takes in sim seqs and names instead of assignments and will run CHMMera directly, instead of CHMMAIRRa
function calculate_plot_four_methods_ROC(sim_seqs::Vector{String}, sim_seq_names::Vector{String}, label, location::String, reference_sets::Dict, reference_set_name::String, prior_probability::Float64, shm1::Float64, shm2::Float64; padding = 0.03, CHMMera_cutoff = 0.95, mutation_probabilities = [0.001, 0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2], vsearch_uchime_cutoff = 0.28, vsearch_uchime_mindiv = 0.8, vsearch_uchime_mindiffs = 3, usearch_uchime2_mindiv = 0., title = "", exclude_methods = String[], receptor = "TCR", gene::Char = 'v', V_seq = "", D_seq = "", J_seq = "", ig_seqtype = "TCR", vsearch_uchime_xn = 8)
    refnames, refseqs = reference_sets[reference_set_name]
    
    # Define all possible methods
    all_methods = ["CHMMAIRRa BW", "CHMMAIRRa DB", "USEARCH uchime2_ref", "VSEARCH uchime_ref"]
    
    # Determine which methods to include
    include_chmmera_bw = !("CHMMAIRRa BW" in exclude_methods)
    include_chmmera_db = !("CHMMAIRRa DB" in exclude_methods)
    include_usearch = !("USEARCH uchime2_ref" in exclude_methods)
    include_vsearch = !("VSEARCH uchime_ref" in exclude_methods)
    
    # Storage for plot data
    plot_fprs = []
    plot_tprs = []
    plot_labels = []
    plot_colors = []
    
    # Compute and store VSEARCH if included
    if include_vsearch
        vsearch_uchime_FPRs, vsearch_uchime_TPRs, vsearch_uchime_scores, vsearch_uchime_cutoffs = vsearch_uchime_ref_ROC_curve(sim_seq_names, degap.(sim_seqs), label, refseqs, mindiv = vsearch_uchime_mindiv, mindiffs = vsearch_uchime_mindiffs, xn = vsearch_uchime_xn)
        push!(plot_fprs, vsearch_uchime_FPRs)
        push!(plot_tprs, vsearch_uchime_TPRs)
        push!(plot_labels, "VSEARCH uchime_ref")
        push!(plot_colors, method2color["VSEARCH uchime_ref"])
    end

    # Compute and store USEARCH if included
    if include_usearch
        usearch_uchime_FPRs, usearch_uchime_TPRs, usearch_uchime_scores, usearch_uchime_cutoffs, usearch_uchime_results = usearch_uchime2_ref_ROC_curve(sim_seq_names, degap.(sim_seqs), label, refseqs, mindiv = usearch_uchime2_mindiv)
        push!(plot_fprs, usearch_uchime_FPRs)
        push!(plot_tprs, usearch_uchime_TPRs)
        push!(plot_labels, "USEARCH uchime2_ref")
        push!(plot_colors, method2color["USEARCH uchime2_ref"])
    end
    

    
    # Compute and store CHMMera BW if included
    if include_chmmera_bw
        CHMMera_BW_FPRs, CHMMera_BW_TPRs, CHMMera_BW_probs, CHMMera_BW_cutoffs = CHMMera_ROC_curve(sim_seqs, label, refseqs, bw = true, prior_probability = prior_probability, mutation_probabilities = mutation_probabilities)
        push!(plot_fprs, CHMMera_BW_FPRs)
        push!(plot_tprs, CHMMera_BW_TPRs)
        push!(plot_labels, "CHMMAIRRa Baum-Welch")
        push!(plot_colors, method2color["CHMMAIRRa BW"])
    end
    
    # Compute and store CHMMera DB if included
    if include_chmmera_db
        CHMMera_DB_FPRs, CHMMera_DB_TPRs, CHMMera_DB_probs, CHMMera_DB_cutoffs = CHMMera_ROC_curve(sim_seqs, label, refseqs, bw = false, prior_probability = prior_probability, mutation_probabilities = mutation_probabilities)
        push!(plot_fprs, CHMMera_DB_FPRs)
        push!(plot_tprs, CHMMera_DB_TPRs)
        push!(plot_labels, "CHMMAIRRa Discretized Bayesian")
        push!(plot_colors, method2color["CHMMAIRRa DB"])
    end
    

    
    # Check that we have at least one method to plot
    if isempty(plot_fprs)
        error("All methods are excluded. At least one method must be included.")
    end
    
    # Create the plot
    p = Plots.plot(plot_fprs, plot_tprs,
        title = title, 
        labels = reshape(plot_labels, 1, :),
        xlabel = "False positive rate", 
        ylabel = "True positive rate", 
        linecolor = reshape(plot_colors, 1, :),
        aspect_ratio = 1.0, 
        markerstrokewidth = 3, 
        legend = :bottomright)
    
    Plots.plot!(p, [0,1], [0,1], color = :black, linestyle = :dash, label = "y = x")
    
    # Storage for annotations
    annotation_x = Float64[]
    annotation_y = Float64[]
    annotation_labels = String[]
    annotation_colors = []
    
    # Add individual cutoff points for CHMMera BW
    if include_chmmera_bw
        cf_ind = findfirst(x->x == CHMMera_cutoff, CHMMera_BW_cutoffs)
        Plots.plot!(p, [CHMMera_BW_FPRs[cf_ind]], [CHMMera_BW_TPRs[cf_ind]], seriestype = :scatter, color = method2color["CHMMAIRRa BW"], label = nothing)
        println("CHMMAIRRa BW cutoff: FPR $(CHMMera_BW_FPRs[cf_ind]) TPR $(CHMMera_BW_TPRs[cf_ind])")
        push!(annotation_x, CHMMera_BW_FPRs[cf_ind])
        push!(annotation_y, CHMMera_BW_TPRs[cf_ind])
        push!(annotation_labels, "P>$(CHMMera_cutoff)")
        push!(annotation_colors, method2color["CHMMAIRRa BW"])
    end
    
    # Add individual cutoff points for CHMMera DB
    if include_chmmera_db
        cs_ind = findfirst(x->x == CHMMera_cutoff, CHMMera_DB_cutoffs)
        Plots.plot!(p, [CHMMera_DB_FPRs[cs_ind]], [CHMMera_DB_TPRs[cs_ind]], seriestype = :scatter, color = method2color["CHMMAIRRa DB"], label = nothing)
        println("CHMMAIRRa DB cutoff: FPR $(CHMMera_DB_FPRs[cs_ind]) TPR $(CHMMera_DB_TPRs[cs_ind])")
        push!(annotation_x, CHMMera_DB_FPRs[cs_ind])
        push!(annotation_y, CHMMera_DB_TPRs[cs_ind])
        push!(annotation_labels, "P>$(CHMMera_cutoff)")
        push!(annotation_colors, method2color["CHMMAIRRa DB"])
    end
    
    # Add USEARCH cutoff points
    if include_usearch
        usearch_cutoff_names = ["sensitive", "balanced", "specific", "high_confidence"]
        for cutoff_name in usearch_cutoff_names
            fpr_val = usearch_uchime_results[cutoff_name].FPR
            tpr_val = usearch_uchime_results[cutoff_name].TPR
            push!(annotation_x, fpr_val)
            push!(annotation_y, tpr_val)
            push!(annotation_labels, cutoff_name)
            push!(annotation_colors, method2color["USEARCH uchime2_ref"])
        end
        
        Plots.plot!(p, [usearch_uchime_results["high_confidence"].FPR, usearch_uchime_results["specific"].FPR, usearch_uchime_results["sensitive"].FPR, usearch_uchime_results["balanced"].FPR],
                [usearch_uchime_results["high_confidence"].TPR, usearch_uchime_results["specific"].TPR, usearch_uchime_results["sensitive"].TPR, usearch_uchime_results["balanced"].TPR], 
                seriestype = :scatter, color = method2color["USEARCH uchime2_ref"], label = nothing)
    end
    
    # Add VSEARCH cutoff points
    if include_vsearch
        u_ind = findfirst(x->x == vsearch_uchime_cutoff, vsearch_uchime_cutoffs)
        if isnothing(u_ind)
            u_ind = length(vsearch_uchime_TPRs)
        end
        Plots.plot!(p, [vsearch_uchime_FPRs[u_ind]], [vsearch_uchime_TPRs[u_ind]], seriestype = :scatter, color = method2color["VSEARCH uchime_ref"], label = nothing)
        println("VSEARCH cutoff: FPR $(vsearch_uchime_FPRs[u_ind]) TPR $(vsearch_uchime_TPRs[u_ind])")
        push!(annotation_x, vsearch_uchime_FPRs[u_ind])
        push!(annotation_y, vsearch_uchime_TPRs[u_ind])
        push!(annotation_labels, "score>$(vsearch_uchime_cutoff)")
        push!(annotation_colors, method2color["VSEARCH uchime_ref"])
    end
    
    # Add annotations with position adjustment
    if !isempty(annotation_x)
        annotation_x_adj = annotation_x .+ padding
        annotation_y_adj = adjust_y_positions(annotation_x_adj, annotation_y, padding = padding)
        for i in 1:length(annotation_x_adj)
            annotate!(p, annotation_x_adj[i], annotation_y_adj[i], Plots.text(annotation_labels[i], annotation_colors[i], :left, 10))
        end
    end
    
    # Print results and AUCs for included methods
    if include_usearch
        print(usearch_uchime_results)
        println("AUC for USEARCH uchime2_ref: $(AUC(usearch_uchime_TPRs,usearch_uchime_FPRs))")
    end
    
    if include_chmmera_bw
        println("AUC for CHMMAIRRa BW: $(AUC(CHMMera_BW_TPRs,CHMMera_BW_FPRs))")
    end
    
    if include_chmmera_db
        println("AUC for CHMMAIRRa DB: $(AUC(CHMMera_DB_TPRs,CHMMera_DB_FPRs))")
    end
    
    if include_vsearch
        println("AUC for VSEARCH uchime_ref: $(AUC(vsearch_uchime_TPRs,vsearch_uchime_FPRs))")
    end
    
    return p
end


# spreads out a set of points based on their y coordinates
# attempts to keep the y order while spreading, so the topmost point remains on top etc
function adjust_y_positions(x::Vector, y::Vector; padding=0.05)
    # adjust from top to bottom to preserve relative order
    y_sortperm = reverse(sortperm(y))
    y = y[y_sortperm]
    positions = [(x[i], y[i]) for i in 1:length(x)]
    adjusted = copy(y)
    overlaps = true
    # keep moving points until none overlap
    while overlaps
        overlaps = false
        for i in 1:length(positions)
            for j in 1:length(positions)
                ydiff = adjusted[i] - adjusted[j]
                # if the points are too close, push the one below down a bit
                # distance to push depends on padding and existing distance between the points
                if (i != j) & (ydiff >=0) & (ydiff < padding)
                    adjusted[j] = adjusted[j] - (rand() * padding * (ydiff + 0.0001) ^ 2)
                    overlaps = true
                end
            end
        end
    end
    unsort_inds = sortperm(y_sortperm)
    return adjusted[unsort_inds]
end



# this one requires IgBLAST assignments as input to run the full CHMMAIRRa method
function calculate_plot_four_methods_ROC(test_sets::DataFrame, location::String, reference_sets::Dict, reference_set_name::String, prior_probability::Float64, shm1::Float64, shm2::Float64; padding = 0.03, CHMMera_cutoff = 0.95, mutation_probabilities = [0.001, 0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2], vsearch_uchime_cutoff = 0.28, title = "", exclude_methods = String[], receptor = "TCR", gene::Char = 'v', V_seq = "", D_seq = "", J_seq = "", ig_seqtype = "TCR", vsearch_uchime_mindiv = 0.8, vsearch_uchime_mindiffs = 3, vsearch_uchime_xn = 8, usearch_uchime2_mindiv = 0.00001)
    curr_simdata = test_sets
    refnames, refseqs = reference_sets[reference_set_name]
    
    # Define all possible methods
    all_methods = ["CHMMAIRRa BW", "CHMMAIRRa DB", "USEARCH uchime2_ref", "VSEARCH uchime_ref"]
    
    # Determine which methods to include
    include_chmmera_bw = !("CHMMAIRRa BW" in exclude_methods)
    include_chmmera_db = !("CHMMAIRRa DB" in exclude_methods)
    include_usearch = !("USEARCH uchime2_ref" in exclude_methods)
    include_vsearch = !("VSEARCH uchime_ref" in exclude_methods)
    
    # Storage for plot data
    plot_fprs = []
    plot_tprs = []
    plot_labels = []
    plot_colors = []


    # Compute and store VSEARCH if included
    if include_vsearch
        vsearch_uchime_FPRs, vsearch_uchime_TPRs, vsearch_uchime_scores, vsearch_uchime_cutoffs = vsearch_uchime_ref_ROC_curve(curr_simdata.sequence_id, degap.(curr_simdata.v_sequence_alignment), curr_simdata.label, refseqs, xn = vsearch_uchime_xn, mindiv = vsearch_uchime_mindiv, mindiffs = vsearch_uchime_mindiffs)
        push!(plot_fprs, vsearch_uchime_FPRs)
        push!(plot_tprs, vsearch_uchime_TPRs)
        push!(plot_labels, "VSEARCH uchime_ref")
        push!(plot_colors, method2color["VSEARCH uchime_ref"])
    end

    # Compute and store USEARCH if included
    if include_usearch
        usearch_uchime_FPRs, usearch_uchime_TPRs, usearch_uchime_scores, usearch_uchime_cutoffs, usearch_uchime_results = usearch_uchime2_ref_ROC_curve(curr_simdata.sequence_id, degap.(curr_simdata.v_sequence_alignment), curr_simdata.label, refseqs, mindiv = usearch_uchime2_mindiv)
        push!(plot_fprs, usearch_uchime_FPRs)
        push!(plot_tprs, usearch_uchime_TPRs)
        push!(plot_labels, "USEARCH uchime2_ref")
        push!(plot_colors, method2color["USEARCH uchime2_ref"])
    end
    
    # Compute and store CHMMera BW if included
    if include_chmmera_bw
        CHMMera_BW_FPRs, CHMMera_BW_TPRs, CHMMera_BW_probs, CHMMera_BW_cutoffs  = CHMMAIRRa_ROC_curve(curr_simdata, curr_simdata.label, refseqs, refnames, receptor = receptor, gene = gene, V_seq = V_seq, D_seq = D_seq, J_seq = J_seq, bw = true, ig_seqtype = ig_seqtype, prior_probability = prior_probability)
        push!(plot_fprs, CHMMera_BW_FPRs)
        push!(plot_tprs, CHMMera_BW_TPRs)
        push!(plot_labels, "CHMMAIRRa Baum-Welch")
        push!(plot_colors, method2color["CHMMAIRRa BW"])
    end
    
    # Compute and store CHMMera DB if included
    if include_chmmera_db
        CHMMera_DB_FPRs, CHMMera_DB_TPRs, CHMMera_DB_probs, CHMMera_DB_cutoffs = CHMMAIRRa_ROC_curve(curr_simdata, curr_simdata.label, refseqs, refnames, receptor = receptor, gene = gene, V_seq = V_seq, D_seq = D_seq, J_seq = J_seq, bw = false, ig_seqtype = ig_seqtype, prior_probability = prior_probability)
        push!(plot_fprs, CHMMera_DB_FPRs)
        push!(plot_tprs, CHMMera_DB_TPRs)
        push!(plot_labels, "CHMMAIRRa Discretized Bayesian")
        push!(plot_colors, method2color["CHMMAIRRa DB"])
    end
    
 
    
    # Check that we have at least one method to plot
    if isempty(plot_fprs)
        error("All methods are excluded. At least one method must be included.")
    end
    
    # Create the plot
    p = Plots.plot(plot_fprs, plot_tprs,
        title = title, 
        labels = reshape(plot_labels, 1, :),
        xlabel = "False positive rate", 
        ylabel = "True positive rate", 
        linecolor = reshape(plot_colors, 1, :),
        aspect_ratio = 1.0, 
        markerstrokewidth = 3, 
        legend = :bottomright)
    
    Plots.plot!(p, [0,1], [0,1], color = :black, linestyle = :dash, label = "y = x")
    
    # Storage for annotations
    annotation_x = Float64[]
    annotation_y = Float64[]
    annotation_labels = String[]
    annotation_colors = []
    
    # Add individual cutoff points for CHMMera BW
    if include_chmmera_bw
        cf_ind = findfirst(x->x == CHMMera_cutoff, CHMMera_BW_cutoffs)
        Plots.plot!(p, [CHMMera_BW_FPRs[cf_ind]], [CHMMera_BW_TPRs[cf_ind]], seriestype = :scatter, color = method2color["CHMMAIRRa BW"], label = nothing)
        println("CHMMAIRRa BW cutoff: FPR $(CHMMera_BW_FPRs[cf_ind]) TPR $(CHMMera_BW_TPRs[cf_ind])")
        push!(annotation_x, CHMMera_BW_FPRs[cf_ind])
        push!(annotation_y, CHMMera_BW_TPRs[cf_ind])
        push!(annotation_labels, "P>$(CHMMera_cutoff)")
        push!(annotation_colors, method2color["CHMMAIRRa BW"])
    end
    
    # Add individual cutoff points for CHMMera DB
    if include_chmmera_db
        cs_ind = findfirst(x->x == CHMMera_cutoff, CHMMera_DB_cutoffs)
        Plots.plot!(p, [CHMMera_DB_FPRs[cs_ind]], [CHMMera_DB_TPRs[cs_ind]], seriestype = :scatter, color = method2color["CHMMAIRRa DB"], label = nothing)
        println("CHMMAIRRa DB cutoff: FPR $(CHMMera_DB_FPRs[cs_ind]) TPR $(CHMMera_DB_TPRs[cs_ind])")
        push!(annotation_x, CHMMera_DB_FPRs[cs_ind])
        push!(annotation_y, CHMMera_DB_TPRs[cs_ind])
        push!(annotation_labels, "P>$(CHMMera_cutoff)")
        push!(annotation_colors, method2color["CHMMAIRRa DB"])
    end
    
    # Add USEARCH cutoff points
    if include_usearch
        usearch_cutoff_names = ["sensitive", "balanced", "specific", "high_confidence"]
        for cutoff_name in usearch_cutoff_names
            fpr_val = usearch_uchime_results[cutoff_name].FPR
            tpr_val = usearch_uchime_results[cutoff_name].TPR
            push!(annotation_x, fpr_val)
            push!(annotation_y, tpr_val)
            push!(annotation_labels, cutoff_name)
            push!(annotation_colors, method2color["USEARCH uchime2_ref"])
        end
        
        Plots.plot!(p, [usearch_uchime_results["high_confidence"].FPR, usearch_uchime_results["specific"].FPR, usearch_uchime_results["sensitive"].FPR, usearch_uchime_results["balanced"].FPR],
                [usearch_uchime_results["high_confidence"].TPR, usearch_uchime_results["specific"].TPR, usearch_uchime_results["sensitive"].TPR, usearch_uchime_results["balanced"].TPR], 
                seriestype = :scatter, color = method2color["USEARCH uchime2_ref"], label = nothing)
    end
    
    # Add VSEARCH cutoff points
    if include_vsearch
        u_ind = findfirst(x->x == vsearch_uchime_cutoff, vsearch_uchime_cutoffs)
        if isnothing(u_ind)
            u_ind = length(vsearch_uchime_TPRs)
        end
        Plots.plot!(p, [vsearch_uchime_FPRs[u_ind]], [vsearch_uchime_TPRs[u_ind]], seriestype = :scatter, color = method2color["VSEARCH uchime_ref"], label = nothing)
        println("VSEARCH cutoff: FPR $(vsearch_uchime_FPRs[u_ind]) TPR $(vsearch_uchime_TPRs[u_ind])")
        push!(annotation_x, vsearch_uchime_FPRs[u_ind])
        push!(annotation_y, vsearch_uchime_TPRs[u_ind])
        push!(annotation_labels, "score>$(vsearch_uchime_cutoff)")
        push!(annotation_colors, method2color["VSEARCH uchime_ref"])
    end
    
    # Add annotations with position adjustment
    if !isempty(annotation_x)
        annotation_x_adj = annotation_x .+ padding
        annotation_y_adj = adjust_y_positions(annotation_x_adj, annotation_y, padding = padding)
        for i in 1:length(annotation_x_adj)
            annotate!(p, annotation_x_adj[i], annotation_y_adj[i], Plots.text(annotation_labels[i], annotation_colors[i], :left, 10))
        end
    end
    
    # Print results and AUCs for included methods
    if include_usearch
        print(usearch_uchime_results)
        println("AUC for USEARCH uchime2_ref: $(AUC(usearch_uchime_TPRs,usearch_uchime_FPRs))")
    end
    
    if include_chmmera_bw
        println("AUC for CHMMAIRRa BW: $(AUC(CHMMera_BW_TPRs,CHMMera_BW_FPRs))")
    end
    
    if include_chmmera_db
        println("AUC for CHMMAIRRa DB: $(AUC(CHMMera_DB_TPRs,CHMMera_DB_FPRs))")
    end
    
    if include_vsearch
        println("AUC for VSEARCH uchime_ref: $(AUC(vsearch_uchime_TPRs,vsearch_uchime_FPRs))")
    end
    
    return p
end

# spreads out a set of points based on their y coordinates
# attempts to keep the y order while spreading, so the topmost point remains on top etc
function adjust_y_positions(x::Vector, y::Vector; padding=0.05)
    # adjust from top to bottom to preserve relative order
    y_sortperm = reverse(sortperm(y))
    y = y[y_sortperm]
    positions = [(x[i], y[i]) for i in 1:length(x)]
    adjusted = copy(y)
    overlaps = true
    # keep moving points until none overlap
    while overlaps
        overlaps = false
        for i in 1:length(positions)
            for j in 1:length(positions)
                ydiff = adjusted[i] - adjusted[j]
                # if the points are too close, push the one below down a bit
                # distance to push depends on padding and existing distance between the points
                if (i != j) & (ydiff >=0) & (ydiff < padding)
                    adjusted[j] = adjusted[j] - (rand() * padding * (ydiff + 0.0001) ^ 2)
                    overlaps = true
                end
            end
        end
    end
    unsort_inds = sortperm(y_sortperm)
    return adjusted[unsort_inds]
end

# plot ROCs for all four methods we're comparing based on the given data, presumably at a specific mutation
function calculate_plot_CHMMAIRRa_ROC(test_sets::DataFrame, location::String, reference_sets::Dict, reference_set_name::String, prior_probability::Float64, shm1::Float64, shm2::Float64; mutation_probabilities = [0.001, 0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2], padding = 0.03, CHMMera_cutoffs = [0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01], plot_sizes = Dict("ROCs" => (500, 500)))
    curr_simdata = test_sets[(test_sets.shm1 .== shm1) .& (test_sets.shm2 .== shm2) .& (test_sets.location .== location),:]
    refnames, refseqs = reference_sets[reference_set_name]
    CHMMera_BW_FPRs, CHMMera_BW_TPRs, CHMMera_BW_probs, CHMMera_BW_cutoffs = CHMMera_ROC_curve(curr_simdata.sequence, curr_simdata.label, refseqs, bw = true, prior_probability = prior_probability, mutation_probabilities = mutation_probabilities)
    CHMMera_DB_FPRs, CHMMera_DB_TPRs, CHMMera_DB_probs, CHMMera_DB_cutoffs = CHMMera_ROC_curve(curr_simdata.sequence, curr_simdata.label, refseqs, bw = false, prior_probability = prior_probability, mutation_probabilities = mutation_probabilities)

    p = Plots.plot([CHMMera_BW_FPRs, CHMMera_DB_FPRs],
    [CHMMera_BW_TPRs, CHMMera_DB_TPRs],
    title = "SHM rate $(Int(floor(shm1 * 100)))%, $(Int(floor(shm2 * 100)))%", labels = ["CHMMAIRRa BW" "CHMMAIRRa DB"],
    xlabel = "False positive rate", ylabel = "True positive rate", linecolor = [method2color["CHMMAIRRa BW"] method2color["CHMMAIRRa DB"]], size = plot_sizes["ROCs"], aspect_ratio = 1.0, markerstrokewidth = 3, legend = :bottomright)


    x = []
    y = []
    labels = []
    colors = []
    Plots.plot!(p, [0,1], [0,1], color = :black, linestyle = :dash, label = "y = x")
    for CHMMera_cutoff in CHMMera_cutoffs
        # add individual cutoff points
        cf_ind = findfirst(x->x == CHMMera_cutoff, CHMMera_BW_cutoffs)
        Plots.plot!(p, [CHMMera_BW_FPRs[cf_ind]], [CHMMera_BW_TPRs[cf_ind]], seriestype = :scatter, color = method2color["CHMMAIRRa BW"], label = nothing)
        println("CHMMAIRRa BW cutoff: FPR $(CHMMera_BW_FPRs[cf_ind]) TPR $(CHMMera_BW_TPRs[cf_ind])")
        # CHMMera cutoff annotation
        cs_ind = findfirst(x->x == CHMMera_cutoff, CHMMera_DB_cutoffs)
        Plots.plot!(p, [CHMMera_DB_FPRs[cs_ind]], [CHMMera_DB_TPRs[cs_ind]], seriestype = :scatter, color = method2color["CHMMAIRRa DB"], label = nothing)
        println("CHMMAIRRa DB cutoff: FPR $(CHMMera_DB_FPRs[cs_ind]) TPR $(CHMMera_DB_TPRs[cs_ind])")
        push!(x, CHMMera_BW_FPRs[cf_ind] + 0.03)
        push!(y, CHMMera_BW_TPRs[cf_ind])
        push!(x, CHMMera_DB_FPRs[cs_ind] + 0.03)
        push!(y, CHMMera_DB_TPRs[cs_ind])
        push!(labels, "P>$(CHMMera_cutoff)")
        push!(labels, "P>$(CHMMera_cutoff)")
        push!(colors, method2color["CHMMAIRRa BW"])
        push!(colors, method2color["CHMMAIRRa DB"])
    end


    y_adj = adjust_y_positions(x, y, padding = padding)
    for i in 1:length(x)
        annotate!(p, x[i], y_adj[i], Plots.text(labels[i], colors[i], :left, 10), )
    end

    println("AUC for CHMMAIRRa BW: $(AUC(CHMMera_BW_TPRs,CHMMera_BW_FPRs))")
    println("AUC for CHMMAIRRa DB: $(AUC(CHMMera_DB_TPRs,CHMMera_DB_FPRs))")
    return p
end

# Benchmark function to compare different shazam implementations
function benchmark_shazam_methods(test_seqs::Vector{String}, shm_rate::Float64; verbose::Bool = true)
    if verbose
        println("Benchmarking shazam methods with $(length(test_seqs)) sequences...")
    end
    
    methods = [
        ("Original (with vector growth)", () -> add_shazam_shm(test_seqs, shm_rate)),
        ("Vectorized (sapply)", () -> add_shazam_shm_vectorized(test_seqs, shm_rate)),
        ("Parallel (4 cores)", () -> add_shazam_shm_parallel(test_seqs, shm_rate, cores=4)),
        ("Uniform random (for comparison)", () -> add_shm(test_seqs, Random.default_rng(), min_shm=shm_rate, max_shm=shm_rate, shm_method="uniform_random"))
    ]
    
    results = []
    for (name, method) in methods
        if verbose
            println("Testing: $name")
        end
        try
            time = @elapsed result = method()
            push!(results, (name, time, length(result)))
            if verbose
                println("  Time: $(round(time, digits=3))s")
            end
        catch e
            if verbose
                println("  Failed: $e")
            end
            push!(results, (name, NaN, 0))
        end
    end
    
    return results
end
