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
    Plots.plot(
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
    Plots.plot!([realParam], seriestype = "vline", margin = 0.3cm, legend = :false)
    ylims!((1E2, 1E5))
end

""" This function calculates the simulated cell number for a pair of drugs in their exact concentration.  """
function calc_cellNum(pDr1, pDr2, g0, controld1, controld2)
    # here we are saying: if the parameters are in their EC50s, what is the Bliss applied to the cell numbers
    # that is the result of each drug on ec50 separately.
    g1d1, g2d1, _ = predict(pDr1, g0, 96.0)
    g1d2, g2d2, _ = predict(pDr2, g0, 96.0)
    normD1 = 1.0 - (g1d1 + g2d1) / (controld1)
    normD2 = 1.0 - (g1d2 + g2d2) / (controld2)
    combin = -(normD1 + normD2 - (normD1 * normD2) - 1.0) * (controld1 + controld2) / 2
    @assert combin >= 0.0

    return combin
end

""" Calculates the difference between the combination and single cell effect. """
function calc_diff(Hillp, Dr1Ind, Dr2Ind, concs, g0, oneortwo)
    effs = getODEparams(Hillp, concs)
    ec501 = EC50_params(Hillp, Dr1Ind)
    ec502 = EC50_params(Hillp, Dr2Ind)
    bliss_comb = Bliss_params_unit(ec501, ec502, hcat(effs[:, 1, Dr1Ind], effs[:, 1, Dr2Ind]))
    g1, g2, _ = predict(bliss_comb, g0, 96.0)
    if oneortwo == 1 # drug 1
        g1_d1, g2_d1, _ = predict(ec501, g0, 96.0)
        return (g1_d1 + g2_d1) - (g1 + g2)
    else # drug 2
        g1_d2, g2_d2, _ = predict(ec502, g0, 96.0)
        return (g1_d2 + g2_d2) - (g1 + g2)
    end
end

import Calculus

function get_derivative(x, Dr1Ind, Dr2Ind, concs, g0, oneortwo)
    fd(x) = calc_diff(x, Dr1Ind, Dr2Ind, concs, g0, oneortwo)
    return Calculus.gradient(fd, x)
end
