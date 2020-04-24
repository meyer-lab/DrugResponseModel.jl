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

    low = [minimum(conc_l), 1e-9, 1e-9, 0.1, 1e-9, 1e-9, 0.0, 0.0, 0.25, 3, 5, 0, 0]
    high = [maximum(conc_l), 3.0, 3.0, 10.0, 3.0, 3.0, 1.0, 1.0, 0.75, 100, 10, 50, 50]

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

""" A function to calculate std and mean of ODE parameters for each drug. """
function mean_std_params(effs1, effs2, effs3)
    meann = ones(9, 8)
    stdd = ones(9, 8)
    for i = 1:8
        for j = 1:9
            meann[j, i] = mean([effs1[j, i], effs2[j, i], effs3[j, i]])
            stdd[j, i] = std([effs1[j, i], effs2[j, i], effs3[j, i]])
        end
    end
    return meann, stdd
end

""" A Function to find mean and std of data in G1 and G2 separately. """
function mean_std_data(G1_1, G1_2, G1_3, G2_1, G2_2, G2_3)

    meanG1 = ones(189, 8)
    meanG2 = ones(189, 8)
    stdG1 = ones(189, 8)
    stdG2 = ones(189, 8)
    for j = 1:8
        for k = 1:189
            meanG1[k, j] = mean([G1_1[k, j], G1_2[k, j], G1_3[k, j]])
            meanG2[k, j] = mean([G2_1[k, j], G2_2[k, j], G2_3[k, j]])
            stdG1[k, j] = std([G1_1[k, j], G1_2[k, j], G1_3[k, j]])
            stdG2[k, j] = std([G2_1[k, j], G2_2[k, j], G2_3[k, j]])
        end
    end
    return meanG1, meanG2, stdG1, stdG2
end

""" A function to predict G1 and G2 for the three replicates. """
function predict_replicates(p1, p2, p3, concs1, g0)
    t = LinRange(0.0, 95.0, 189)
    G1_1 = ones(189, 8)
    G2_1 = ones(189, 8)
    G1_2 = ones(189, 8)
    G2_2 = ones(189, 8)
    G1_3 = ones(189, 8)
    G2_3 = ones(189, 8)

    for i = 1:8 # concentration number
        G1_1[:, i], G2_1[:, i], _ = predict(p1[:, i], g0, t)
        G1_2[:, i], G2_2[:, i], _ = predict(p2[:, i], g0, t)
        G1_3[:, i], G2_3[:, i], _ = predict(p3[:, i], g0, t)
    end

    return G1_1, G2_1, G1_2, G2_2, G1_3, G2_3 # all simulation
end

""" A function to plot one of the concentrations for the three replicates with ribbon of std data. """
function plot_reps_ribbon(G1_1, G1_2, G1_3, G2_1, G2_2, G2_3, meang1, meang2, stdg1, stdg2, conc, legend)
    time = LinRange(0.0, 95.0, 189)
    title = string(conc, " nM")
    plot(time, meang1; ribbon = stdg1, color = 6, label = "", xlabel = "time [hr]", ylabel = "cell number", alpha = 0.05, legend = legend)
    plot!(time, G1_1, label = "G1", color = 6)
    plot!(time, G1_2, label = "", color = 6)
    plot!(time, G1_3, label = "", color = 6)
    plot!(time, meang2; ribbon = stdg2, color = 7, label = "", alpha = 0.05)
    plot!(time, G2_1, label = "G2", color = 7)
    plot!(time, G2_2, label = "", color = 7)
    plot!(time, G2_3, label = "", color = 7)
    plot!(annotation = [(45, 40, text(title, 8))])
    ylims!((0.0, 45))
end
