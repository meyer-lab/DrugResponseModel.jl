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
    @assert length(params) == 13
    b = copy(params)
    
    convRange = 10 .^ range(-1, stop = 1, length = 101)
    results = zeros(length(convRange), length(params))
    paramRanges = zeros(length(convRange), length(params))

    for k = 1:length(params)
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
