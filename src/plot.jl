gr()
""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

#-------------------------- plot for the G1 G2 correlation ---------------------------#
function correlationPlot(g1::Matrix, g2::Matrix, labels::Array, xlabel::String, ylabel::String, ymax::Int)
    theme(:mute)
    pl = [
        scatter(
            g1[:, i],
            g2[:, i],
            title = labels[i],
            xlabel = xlabel,
            ylabel = ylabel,
            alpha = 0.6,
#             color = [:gray],
            line = :dot,
            marker = (:dot, 1.5),
        )
        for i = 1:8
    ]
    plot(pl..., legend = :false, layout = (2, 4), fmt = :png)
    plot!(size = (1200, 600), margin = 0.4cm, dpi = 150)
    ylims!((0, ymax))
    xlims!((0, ymax))
end


function plot_parameters(conc_l, parameters, costResults, paramRange)

    sd = errorbar_params(conc_l, costResults, paramRange)
    
    conc = log.(conc_l)
    theme(:mute)
    p1 = plot(
        conc,
        parameters[1, :],
        xlabel = "log drug conc. [nM]",
        label = "",
        lw = 2.0,
        alpha = 0.6,
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "alpha",
        ribbon = (sd[1,:])
    )
    ylims!(0.0, 1.2 * maximum(parameters[1, :]))

    p2 = plot(
        conc,
        parameters[2, :],
        xlabel = "log drug conc. [nM]",
        label = "",
        lw = 2.0,
        alpha = 0.6,
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "beta",
        ribbon = (sd[2,:])
    )
    ylims!(0.0, 1.2 * maximum(parameters[2, :]))

    maxDeath = maximum(parameters[3:4, :])

    p3 = plot(
        conc,
        parameters[3, :],
        xlabel = "log drug conc. [nM]",
        label = "",
        lw = 2.0,
        alpha = 0.6,
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "gamma1",
        ribbon = (sd[3,:])
    )
    ylims!(0.0, 1.2 * maxDeath)

    p4 = plot(
        conc,
        parameters[4, :],
        xlabel = "log drug conc. [nM]",
        label = "",
        lw = 2.0,
        alpha = 0.6,
#         color = [:black :gray],
        line = (:dot, 1),
        marker = ([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "gamma2",
        ribbon = (sd[4,:])
    )
    ylims!(0.0, 1.2 * maxDeath)

    plot(p1, p2, p3, p4)
    plot!(size = (600, 400), margin = 0.4cm, dpi = 150)
end

""" Use the results from sensitivity analysis to plot a range for parameters vs. concentrations. """
function errorbar_params(concentrations, results, paramRange)
    newparam = zeros(3, 11)
    for j in 1:11
        ind_min = argmin(results[:, j])
        if (ind_min <= 5 || ind_min >=95)
            ind_min = 50
        end
        newparam[:, j] = paramRange[(ind_min - 5):5:(ind_min + 5), j]
    end
    odeparams = zeros(7,8,3)
    for k in 1:3
        odeparams[:, :, k] = getODEparams(newparam[k, :], concentrations)
    end
    sd = odeparams[:, :, 3] .- odeparams[:, :, 1]
    return sd
end
