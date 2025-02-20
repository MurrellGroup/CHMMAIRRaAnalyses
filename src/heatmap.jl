using CSV, DataFrames, CairoMakie, Plots, Clustering

"""
Calculate the actual co occurrence frequency matrix from a dataframe of chimera V gene co occurrences
Counts are normalized per library, then averaged across libraries
Matrix is ordered using the genotype variable
"""
function co_occurrence_mat(df, genotype)
    name2ind = Dict([(genotype[i], i) for i in 1:length(genotype)])
    order_known = collect(keys(name2ind))

    case_freqmats = Dict()
    for case_df in groupby(df, :case)
        mat = zeros(length(genotype), length(genotype))
        recombinatin_arrays = parse_recombinations_string.(case_df.recombinations_degapped)
        for recombination_array in recombinatin_arrays
            for recombination in recombination_array
                left_gene, right_gene = split(recombination.left_allele, "*")[1], split(recombination.right_allele, "*")[1]
                if haskey(name2ind, left_gene) && haskey(name2ind, right_gene)
                    mat[name2ind[left_gene], name2ind[right_gene]] += 1
                end
            end
        end
        case_freqmats[case_df.case[1]] = mat
    end
    # average recombination frequencies across donors
    mat = zeros(length(genotype), length(genotype))
    for name1 in genotype
        for name2 in genotype
            mat[name2ind[name1], name2ind[name2]] = mean([case_freqmats[case][name2ind[name1], name2ind[name2]] for case in keys(case_freqmats)])
        end
    end
    mat = mat ./ sum(mat)
    return mat
end

"""
Calculate the expected frequency of chimerism from one V gene to another V gene using repertoire only V gene frequencies
Matrix is ordered using the genotype variable
"""
function null_mat(df, col, genotype)
    # frequency per case
    totals = combine(groupby(df, [col, "case"]), col => length => :n)
    case_cts = combine(groupby(df, ["case"]), :case => length => :ncase)
    totals = leftjoin(totals, case_cts, on = [:case])
    totals[!,"frequency"] .= totals[!,"n"] ./ sum(totals[!,"ncase"])

    name2ind = Dict([(genotype[i], i) for i in 1:length(genotype)])
    order_known = collect(keys(name2ind))
    calls = unique(df[!,col])
    @warn "$(setdiff(calls, order_known)) have no known order and will be ignored"
    df = df[map(x-> (x[col] in order_known) , eachrow(df)),:]
    # multiply frequencies per case
    case_freqmats = Dict()
    for case_totals in groupby(totals, :case)
        name2frequency = Dict([(row[col], row.frequency) for row in eachrow(case_totals)])
        mat = zeros(length(genotype), length(genotype))
        for name1 in genotype
            for name2 in genotype
                mat[name2ind[name1], name2ind[name2]] = get(name2frequency, name1, 0) * get(name2frequency, name2, 0)
            end
        end
        case_freqmats[case_totals.case[1]] = mat
    end
    # average recombination frequencies across donors
    mat = zeros(length(genotype), length(genotype))
    for name1 in genotype
        for name2 in genotype
            mat[name2ind[name1], name2ind[name2]] = mean([case_freqmats[case][name2ind[name1], name2ind[name2]] for case in keys(case_freqmats)])
        end
    end
    mat = mat ./ sum(mat)
    return mat
end

"""
Plot co occurrence frequency heatmap
Note the function takes the root of the frequencies
"""
function cocc_heatmap(cocc, order; root = 2, plot_size = (900,900), m = nothing, title = "", titlefontsize = 5, transform_f = x->x^(1/2))
    cocc = transform_f.(cocc)
    m = m == nothing ? maximum(cocc) : m
    p = Plots.plot(cocc, seriestype = :heatmap, c=cgrad([:white, :blue]),
            xticks = (1:1:length(order), order),
            yticks = (1:1:length(order), order),
            xtickfontsize = 5,
            ytickfontsize = 5,
            xrotation = 45,
            size = plot_size,
            tick_direction = :none,
            clims=(0, m),
            bottom_margin=10mm,
            title = title,
            titlefontsize=titlefontsize)
    m, n = size(cocc)
    vline!(p, 0.5:(n+0.5), c=:lightgrey, label = false)
    hline!(p, 0.5:(m+0.5), c=:lightgrey, label = false)
end

# I had difficulty making an independent colorbar, so I created a dumb work around
function manual_vertical_colorbar(min_::Float64, max_::Float64; titlefontsize::Int64 = 5)
	digits = max_ > 1 ? 1 : 2
	i = collect(min_:.0001:max_)
	m = zeros(length(i), 1)
	m[:,1] .= i
	p = Plots.plot(m, seriestype = :heatmap, clims = (min_, max_),
				yticks = false, legend = false, c=cgrad([:white, :blue]),
				yaxis = false, xaxis = false, top_margin = 0mm, bottom_margin = 1mm, right_margin = 3mm,
                title = "√freq", titlefontsize = titlefontsize - 1)
	annotate!(p, 1.7, 0, text(round(min_, digits = digits), :black, :left, 5))
	annotate!(p, 1.7, length(i) / 2,  text(round(.5 * max_, digits = digits), :black, :left, 5))
	annotate!(p, 1.7, length(i),  text(round(max_, digits = digits), :black, :left, 5))
	return p
end

# from https://github.com/MakieOrg/Makie.jl/issues/398
function treepositions(hc::Hclust; useheight = true, orientation = :vertical)
    order = StatsBase.indexmap(hc.order)
    nodepos = Dict(-i => (float(order[i]), 0.0) for i in hc.order)
    xs = []
    ys = []
    for i in 1:size(hc.merges, 1)
        x1, y1 = nodepos[hc.merges[i, 1]]
        x2, y2 = nodepos[hc.merges[i, 2]]
        xpos = (x1 + x2) / 2
        ypos = useheight ?  hc.heights[i] : (max(y1, y2) + 1)
        nodepos[i] = (xpos, ypos)
        push!(xs, [x1, x1, x2, x2])
        push!(ys, [y1, ypos, ypos, y2])
    end
    if orientation == :horizontal
        return ys, xs
    else
        return xs, ys
    end
end

function cocc_heatmap_dendograms_CairoMakie(cocc, hc, order; plot_size = (500, 500), root = 2, title = "", titlefontsize = 5, transform_f = x->x^(1/2), heatmapticklabelsize = 20, dendogramticklabelsize = 15, colorbarticklabelsize = 25, colorbarlabelsize = 20, heatmapaxislabelsize = 20)

    cocc = transform_f.(cocc)
    n, m = size(cocc)

    f = CairoMakie.Figure(size = plot_size, figure_padding = (30,70,30,30))
    Label(f[1,:], text = title, fontsize = titlefontsize, valign = :bottom)

    # need to enforce dendogram limits to align dendogram with heatmap
    dend_x_coords = collect(Iterators.flatten([x for (x,y) in zip(treepositions(hc)...)]))
    dend_x_max, dend_x_min = maximum(dend_x_coords) + 0.5, minimum(dend_x_coords) - 0.5

    top_dendogram = CairoMakie.Axis(f[2, 1], xticksize = 0, yticklabelsize = dendogramticklabelsize, limits = (dend_x_min, dend_x_max, nothing, nothing), xgridvisible = false, ygridvisible = false, xticklabelsize = 0)
    heatmap = CairoMakie.Axis(f[3, 1], xticks = (1:m, order), yticks = (1:m, order), xticklabelrotation = π / 2, xticklabelsize = heatmapticklabelsize, yticklabelsize = heatmapticklabelsize, xlabel = "Right parent gene", ylabel = "Left parent gene", xlabelsize = heatmapaxislabelsize, ylabelsize = heatmapaxislabelsize)
    right_dendogram = CairoMakie.Axis(f[3, 2], yticksvisible = false, xticklabelrotation = π / 2, xticklabelsize = dendogramticklabelsize, limits = (nothing, nothing, dend_x_min - m, dend_x_max - m), xgridvisible = false, ygridvisible = false, yticklabelsize = 0)

    CairoMakie.heatmap!(heatmap, cocc, colormap = cgrad([:white, :blue]))
    # diagonal line
    CairoMakie.lines!(heatmap, [0.5, m + 0.5], [0.5, m + 0.5], color = :grey)
    # grid lines
    CairoMakie.hlines!(heatmap, .5:1:m, color = :grey)
    CairoMakie.vlines!(heatmap, .5:1:n, color = :grey)

    # dendograms
    for (x, y) in zip(treepositions(hc)...)
        CairoMakie.lines!(top_dendogram, x, y; color = :black)
        CairoMakie.lines!(right_dendogram, y, - m .+ x; color = :black)
    end
    hidespines!(top_dendogram)
    hidespines!(right_dendogram)

    # colorbar for the heatmap
    cb = Colorbar(f[3,3], limits = (0, maximum(cocc)), colormap = cgrad([:white, :blue]), halign = :left, size = 15, ticklabelsize = colorbarticklabelsize, labelsize = colorbarlabelsize, label = "√Frequency")
    cb.alignmode = Mixed(right = 0)

    # remove the gap between the heatmap and the dendograms
    colgap!(f.layout, 1, 3)
    colgap!(f.layout, 2, 3)

    rowgap!(f.layout, 1, 8)
    rowgap!(f.layout, 2, 3)

    # set relative size of heatmap
    rowsize!(f.layout, 1, Relative(.1))
    colsize!(f.layout, 1, Relative(.8))
    rowsize!(f.layout, 2, Relative(.1))
    colsize!(f.layout, 2, Relative(.1))
    rowsize!(f.layout, 3, Relative(.8))
    colsize!(f.layout, 3, Relative(.1))

    # make the heatmap square
    rowsize!(f.layout, 3, Aspect(1, 1))

    current_figure()
end