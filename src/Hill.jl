"""
This file fits Hill function to the parameters
"""

""" This functions takes in hill parameters for all the concentrations and calculates
DDE parameters, passes them to residual function and based off of these, optimizes the model
and estimates hill parameters. """
function residHill(hillParams::Vector, concentrations::Vector, g1::Matrix, g2::Matrix)
    res = 0.0
    params = getODEparams(hillParams, concentrations)

    # Solve each concentration separately
    for ii = 1:length(concentrations)
        resTemp = cost(params[1:9, ii], g1[:, ii], g2[:, ii])

        res += resTemp
    end

    return res
end

""" Hill optimization function. """
function optimize_hill(conc_l::Vector, g1::Matrix, g2::Matrix; maxstep = 1E5)
    hillCost(hillParams) = residHill(hillParams, conc_l, g1, g2)

    low = [minimum(conc_l), 1e-9, 1e-9, 0.1, 1e-9, 1e-9, 0.0, 0.0, 0.45, 2, 10, 0, 0]
    high = [maximum(conc_l), 3.0, 3.0, 10.0, 3.0, 3.0, 1.0, 1.0, 0.55, 60, 180, 50, 50]

    results_ode = bboptimize(
        hillCost;
        SearchRange = collect(zip(low, high)),
        NumDimensions = length(low),
        TraceMode = :verbose,
        TraceInterval = 100,
        MaxSteps = maxstep,
    )

    return best_fitness(results_ode), best_candidate(results_ode)
end

""" A function to convert the estimated hill parameters back to ODE parameters. """
function getODEparams(p::Vector, concentrations::Vector{Float64})
    effects = Matrix{eltype(p)}(undef, 9, length(concentrations))

    # Scaled drug effect
    xx = 1.0 ./ (1.0 .+ (p[1] ./ (concentrations .+ eps())) .^ p[4])

    # [EC50, left, right, steepness]
    effects[1, :] = p[2] .+ (p[3] - p[2]) .* xx
    effects[2, :] = p[5] .+ (p[6] - p[5]) .* xx
    effects[3, :] = p[7] .* xx
    effects[4, :] = p[8] .* xx
    effects[5, :] .= p[9]
    effects[6, :] .= floor(p[10])
    effects[7, :] .= floor(p[11])
    effects[8, :] .= floor(p[12])
    effects[9, :] .= floor(p[13])

    return effects
end

""" To find the sensitivity of the model to a parameter. """
function sensitivity(params::Vector, paramRange::Vector, conc::Vector, i::Int, g1::Matrix, g2::Matrix)
    result = zeros(length(paramRange))
    for j = 1:length(paramRange)
        temp = copy(params)
        temp[i] = paramRange[j]
        result[j] = residHill(temp, conc, g1, g2)
    end
    return result
end

""" Calculate the sensitivity to all parameters. """
function allSensitivity(params::Vector, conc_l::Vector, g1::Matrix, g2::Matrix)
    b = copy(params)
    convRange = 10 .^ (range(-1, stop = 1, length = 101))
    results = zeros(length(convRange), 11)
    paramRanges = zeros(length(convRange), 11)

    for k = 1:11
        paramRanges[:, k] = b[k] .* convRange
        results[:, k] = sensitivity(b, paramRanges[:, k], conc_l, k, g1, g2)
    end

    return results, paramRanges
end

""" Plots the sensitivity for a parameter with a vertical line of the real value of the parameter."""
function plotUnitSensitivity(paramRange, result, realParam, i)
    label = [
        "EC50",
        "min alpha",
        "max alpha",
        "steepnsess",
        "min beta",
        "max beta",
        "max gamma G1",
        "max gamma G2",
        "% in G1",
        "# of G1 species",
        "# of G2 species",
    ]
    plot(
        paramRange,
        result,
        legend = :false,
        xlabel = string("[log] ", label[i], " range"),
        ylabel = "[log] cost",
        yscale = :log10,
        xscale = :log10,
        xaxis = :log10,
        yaxis = :log10,
    )
    plot!([realParam], seriestype = "vline", margin = 0.3cm, legend = :false)
    ylims!((1E2, 1E4))
end


""" Calculate the # of cells in G1 for a set of parameters and T """
function numcells(params, g0, T)
    @assert(all(params .>= 0.0), "negative params $params")
    t = LinRange(0.0, 95.5, 192)
    G1, G2 = predict(params, g0, t)

    @assert(all(G1[2:end] .>= 0.0), "negative cell number in G1 $G1")
    @assert(all(G2[2:end] .>= 0.0), "negative cell number in G2 $G2")
    return G1[T] + G2[T]
end
