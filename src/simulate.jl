using CSV, DataFrames, Random, StatsBase

# samples n alleles per gene from the reference set
function simulate_genotype(refnames, refseqs, rng; n_alleles_per_gene = 1)
    genenames = map(refname -> split(refname, "*")[1], refnames)
    unique_genenames = unique(genenames)
    u_gene_inds = vcat([sample(rng, findall(genenames .== genename), n_alleles_per_gene) for genename in unique_genenames]...)
    # take one allele per gene
    refnames, refseqs = refnames[u_gene_inds], refseqs[u_gene_inds]
    refnames, refseqs = mafft_wrapper(refseqs, refnames)
    return refnames, refseqs
end

function random_seqpair(names::Vector{String}, seqs::Vector{String}, rng)
    rand_inds = sample(rng, 1:length(names), 2, replace = false)
    return names[rand_inds], seqs[rand_inds]
end

function get_breakpoint(seq::String, rng; min_pos::Union{Int, Float64} = 0.0, max_pos::Union{Int, Float64} = 1.0)
    if ismissing(max_pos)
        max_pos = length(seq)
    end
    if (typeof(min_pos) == Float64) & (min_pos < 1)
        min_pos = floor(min_pos * length(seq))
    end
    if (typeof(max_pos) == Float64) & (max_pos <= 1)
        max_pos = floor(max_pos * length(seq))
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

# pick n random sequences from the reference set and add shm to them
# shm rate chosen uniformly between min_shm and max_shm
# assumes sequences are already aligned
function random_nonchimeras(names::Vector{String}, seqs::Vector{String}, rng; n::Int = 1, min_shm::Float64 = 0.2, max_shm::Float64 = 0.2)
    if length(unique(length.(seqs))) != 1
        error("All sequences must be the same length. Please align them and try again.")
    end
    rand_inds = sample(rng, 1:length(names), n, replace = true)
    rand_names, rand_seqs = names[rand_inds], seqs[rand_inds]
    rand_seqs = map(rand_seq -> add_random_shm(rand_seq, rng, min_shm = min_shm, max_shm = max_shm), rand_seqs) # add shm
    rand_names = ["$(rand_names[i])_$(i)" for i in 1:length(rand_names)]
    return rand_names, rand_seqs
end

# Generate n chimeras by taking random sequences, cutting at a random breakpoint between min_pos and max_pos, and stitching together the cut sequences
# assumes sequences are already aligned
function random_chimeras(names::Vector{String}, seqs::Vector{String}, rng; min_pos::Union{Int, Float64} = 0.0, max_pos::Union{Int, Float64} = 1.0, n::Int = 1, min_shm1::Float64 = 0.0, max_shm1::Float64 = 0.2, min_shm2::Float64 = 0.0, max_shm2::Float64 = 0.2)
    if length(unique(length.(seqs))) != 1
        error("All sequences must be the same length. Please align them and try again.")
    end
    rand_refs_list = [random_seqpair(names, seqs, rng) for i in 1:n] # choose seqeunces
    rand_refs_list = map(rand_refs -> [rand_refs[1], [add_random_shm(rand_refs[2][1], rng, min_shm = min_shm1, max_shm = max_shm1), add_random_shm(rand_refs[2][2], rng, min_shm = min_shm2, max_shm = max_shm2)]], rand_refs_list) # add shm
    breakpoint_positions = [get_breakpoint(rand_refs[2][1], rng, min_pos = min_pos, max_pos = max_pos) for rand_refs in rand_refs_list] # choose breakpoints
    chimeric_seqs = [string(rand_refs[2][1][1:curr_breakpoint_position], rand_refs[2][2][curr_breakpoint_position + 1:end]) for (rand_refs, curr_breakpoint_position) in zip(rand_refs_list, breakpoint_positions)] # concatenate seq fragments
    chimeric_names = ["$(rand_refs[1][1])_1to$(curr_breakpoint_position)_$(rand_refs[1][2])_$(curr_breakpoint_position)toEnd_$(i)" for (i, rand_refs, curr_breakpoint_position) in zip(1:length(rand_refs_list), rand_refs_list, breakpoint_positions)]
    return chimeric_names, chimeric_seqs, [el[1] for el in rand_refs_list], [el[2] for el in rand_refs_list], breakpoint_positions
end

function usearch_uchime2_ref_wrapper(query_names::Vector{String}, query_seqs::Vector{String}, db_seqs::Vector{String}; search::String = "usearch", mode = "sensitive")
    """
    Runs uchime_ref and returns sequence scores
    Note that query_names must be unique because we need them to sort the uchime_ref output
    """
    local res
    mktempdir() do dir
        # run uchime_ref
        db_fasta_path = joinpath(dir, "V.fasta")
        query_fasta_path = joinpath(dir, "queries.fasta")
        outpath = joinpath(dir, "out.txt")
        write_fasta(db_fasta_path, degap.(db_seqs))
        write_fasta(query_fasta_path, degap.(query_seqs), seq_names = query_names)
        cmd = `usearch -uchime2_ref $(query_fasta_path) -db $(db_fasta_path) -uchimeout $(outpath) -strand plus -mode $(mode)`
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

function parse_attributes(attributes::String)
    attribute_dict::Dict{String, Any} = Dict("dqt" => missing, "dqm" => missing, "div" => missing, "L_Y" => missing, "L_N" => missing, "L_A" => missing, "R_Y" => missing, "R_N" => missing, "R_A" => missing, "L" => missing, "R" => missing, "why" => missing)
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
    return DataFrame(attribute_dict)
end

function vsearch_uchime_ref_wrapper(query_names::Vector{String}, query_seqs::Vector{String}, db_seqs::Vector{String}; vsearch::String = "vsearch")
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
        cmd = `conda run -n vsearch vsearch --uchime_ref $(query_fasta_path) --uchimeout $(outpath) --fasta_score --db $(db_fasta_path)`
        println(cmd)
        run(cmd)
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

function usearch_uchime2_ref_ROC_curve(query_names::Vector{String}, query_seqs::Vector{String}, labels, ref_seqs::Vector{String}; mindiffs = 2, mindivt = 0.5, xa = 1, xn = 4, modes = ["high_confidence", "specific", "balanced", "sensitive"])
    """
    Generate an ROC curve plot for a set of query sequences using the CHMMera method by varying the chimeric probability threshold
    """
    calls = Dict()
    # we'll want the last output from the loop, so we need to initialize the variables outside of the loop
    TPRs, FPRs, scores, cutoffs = zeros(0), zeros(0), zeros(0), zeros(0)
    for mode in modes
        res = usearch_uchime2_ref_wrapper(query_names, query_seqs, ref_seqs, mode = mode)
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

function vsearch_uchime_ref_ROC_curve(query_names::Vector{String}, query_seqs::Vector{String}, labels, ref_seqs::Vector{String}; cutoff_interval = 0.01)
    """
    Generate an ROC curve plot for a set of query sequences using the CHMMera method by varying the chimeric probability threshold
    """
    # get uchime_ref scores
    scores = vsearch_uchime_ref_wrapper(query_names, query_seqs, ref_seqs)
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

function AUC(TPRs::Vector{Float64}, FPRs::Vector{Float64})
    """Calculate AUC using trapezoids"""
    AUC = 0
    for i in 2:length(TPRs)
        AUC = AUC +  ((FPRs[i] - FPRs[i - 1]) * ((TPRs[i] + TPRs[i - 1]) / 2))
    end
    return abs(AUC)
end



# plot ROCs for all four methods we're comparing based on the given data, presumably at a specific mutation
function calculate_plot_four_methods_ROC(test_sets::DataFrame, location::String, reference_sets::Dict, reference_set_name::String, prior_probability::Float64, shm1::Float64, shm2::Float64; padding = 0.03, CHMMera_cutoff = 0.95, mutation_probabilities = [0.001, 0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2], vsearch_uchime_cutoff = 0.28, title = "")
    curr_simdata = test_sets[(test_sets.shm1 .== shm1) .& (test_sets.shm2 .== shm2) .& (test_sets.location .== location) .& (test_sets.refset_name .== reference_set_name),:]
    refnames, refseqs = reference_sets[reference_set_name]
    CHMMera_BW_FPRs, CHMMera_BW_TPRs, CHMMera_BW_probs, CHMMera_BW_cutoffs = CHMMera_ROC_curve(curr_simdata.sequence, curr_simdata.label, refseqs, bw = true, prior_probability = prior_probability, mutation_probabilities = mutation_probabilities)
    CHMMera_DB_FPRs, CHMMera_DB_TPRs, CHMMera_DB_probs, CHMMera_DB_cutoffs = CHMMera_ROC_curve(curr_simdata.sequence, curr_simdata.label, refseqs, bw = false, prior_probability = prior_probability, mutation_probabilities = mutation_probabilities)
    usearch_uchime_FPRs, usearch_uchime_TPRs, usearch_uchime_scores, usearch_uchime_cutoffs, usearch_uchime_results = usearch_uchime2_ref_ROC_curve(curr_simdata.sequence_id, curr_simdata.sequence, curr_simdata.label, refseqs)
    vsearch_uchime_FPRs, vsearch_uchime_TPRs, vsearch_uchime_scores, vsearch_uchime_cutoffs = vsearch_uchime_ref_ROC_curve(curr_simdata.sequence_id, curr_simdata.sequence, curr_simdata.label, refseqs)

    p = Plots.plot([CHMMera_BW_FPRs, CHMMera_DB_FPRs, usearch_uchime_FPRs, vsearch_uchime_FPRs],
    [CHMMera_BW_TPRs, CHMMera_DB_TPRs, usearch_uchime_TPRs, vsearch_uchime_TPRs],
    title = title, labels = ["CHMMAIRRa BW" "CHMMAIRRa DB" "USEARCH uchime2_ref" "VSEARCH uchime_ref"],
    xlabel = "False positive rate", ylabel = "True positive rate", linecolor = [method2color["CHMMAIRRa BW"] method2color["CHMMAIRRa DB"] method2color["USEARCH uchime2_ref"] method2color["VSEARCH uchime_ref"]], aspect_ratio = 1.0, markerstrokewidth = 3, legend = :bottomright)


    Plots.plot!(p, [0,1], [0,1], color = :black, linestyle = :dash, label = "y = x")
    # add individual cutoff points
    cf_ind = findfirst(x->x == CHMMera_cutoff, CHMMera_BW_cutoffs)
    Plots.plot!(p, [CHMMera_BW_FPRs[cf_ind]], [CHMMera_BW_TPRs[cf_ind]], seriestype = :scatter, color = method2color["CHMMAIRRa BW"], label = nothing)
    println("CHMMAIRRa BW cutoff: FPR $(CHMMera_BW_FPRs[cf_ind]) TPR $(CHMMera_BW_TPRs[cf_ind])")
    # CHMMera cutoff annotation
    cs_ind = findfirst(x->x == CHMMera_cutoff, CHMMera_DB_cutoffs)
    Plots.plot!(p, [CHMMera_DB_FPRs[cs_ind]], [CHMMera_DB_TPRs[cs_ind]], seriestype = :scatter, color = method2color["CHMMAIRRa DB"], label = nothing)
    println("CHMMAIRRa DB cutoff: FPR $(CHMMera_DB_FPRs[cs_ind]) TPR $(CHMMera_DB_TPRs[cs_ind])")

    # USEARCH uchime2_ref cutoff annotation
    Plots.plot!(p, [usearch_uchime_results["high_confidence"].FPR, usearch_uchime_results["specific"].FPR, usearch_uchime_results["sensitive"].FPR, usearch_uchime_results["balanced"].FPR],
            [usearch_uchime_results["high_confidence"].TPR, usearch_uchime_results["specific"].TPR, usearch_uchime_results["sensitive"].TPR, usearch_uchime_results["balanced"].TPR], seriestype = :scatter, color = method2color["USEARCH uchime2_ref"], label = nothing)

    # VSEARCH uchime_ref cutoff annotation
    u_ind = findfirst(x->x == vsearch_uchime_cutoff, vsearch_uchime_cutoffs)
    if isnothing(u_ind)
        u_ind = length(vsearch_uchime_TPRs)
    end
    Plots.plot!(p, [vsearch_uchime_FPRs[u_ind]], [vsearch_uchime_TPRs[u_ind]], seriestype = :scatter, color = method2color["VSEARCH uchime_ref"], label = nothing)
    println("VSEARCH cutoff: FPR $(vsearch_uchime_FPRs[u_ind]) TPR $(vsearch_uchime_TPRs[u_ind])")

    x = [CHMMera_BW_FPRs[cf_ind], CHMMera_DB_FPRs[cs_ind], usearch_uchime_results["sensitive"].FPR, usearch_uchime_results["balanced"].FPR, usearch_uchime_results["specific"].FPR, usearch_uchime_results["high_confidence"].FPR, vsearch_uchime_FPRs[u_ind]] .+ 0.03
    y = [CHMMera_BW_TPRs[cf_ind], CHMMera_DB_TPRs[cs_ind], usearch_uchime_results["sensitive"].TPR, usearch_uchime_results["balanced"].TPR, usearch_uchime_results["specific"].TPR, usearch_uchime_results["high_confidence"].TPR, vsearch_uchime_TPRs[u_ind]]
    labels = ["P>$(CHMMera_cutoff)", "P>$(CHMMera_cutoff)", "sensitive", "balanced", "specific", "high_confidence", "score>$(vsearch_uchime_cutoff)"]
    colors = [method2color["CHMMAIRRa BW"], method2color["CHMMAIRRa DB"], method2color["USEARCH uchime2_ref"], method2color["USEARCH uchime2_ref"], method2color["USEARCH uchime2_ref"], method2color["USEARCH uchime2_ref"], method2color["VSEARCH uchime_ref"]]

    y_adj = adjust_y_positions(x, y, padding = padding)
    for i in 1:length(x)
        annotate!(p, x[i], y_adj[i], text(labels[i], colors[i], :left, 10), )
    end

    print(usearch_uchime_results)

    println("AUC for CHMMAIRRa BW: $(AUC(CHMMera_BW_TPRs,CHMMera_BW_FPRs))")
    println("AUC for CHMMAIRRa DB: $(AUC(CHMMera_DB_TPRs,CHMMera_DB_FPRs))")
    println("AUC for USEARCH uchime2_ref: $(AUC(usearch_uchime_TPRs,usearch_uchime_FPRs))")
    println("AUC for VSEARCH uchime_ref: $(AUC(vsearch_uchime_TPRs,vsearch_uchime_FPRs))")
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
        annotate!(p, x[i], y_adj[i], text(labels[i], colors[i], :left, 10), )
    end

    println("AUC for CHMMAIRRa BW: $(AUC(CHMMera_BW_TPRs,CHMMera_BW_FPRs))")
    println("AUC for CHMMAIRRa DB: $(AUC(CHMMera_DB_TPRs,CHMMera_DB_FPRs))")
    return p
end
