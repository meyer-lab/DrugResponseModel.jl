gr()
""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

#-------------------------- plot for the G1 G2 correlation ---------------------------#
function correlationPlot(g1::Matrix, g2::Matrix, labels::Array, xlabel::String, ylabel::String, ymax::Int)
    pl = [scatter(g1[:, i], g2[:, i], title=labels[i], xlabel=xlabel, ylabel=ylabel, alpha=0.6, color=[:gray], line=:dot, marker=(:dot, 1.5)) for i in 1:8]
    plot(pl..., legend=:false, layout=(2,4), fmt = :png)
    plot!(size=(1200,600), dpi=150)
    ylims!((0, ymax))
    xlims!((0, ymax))
end


function plot_parameters(conc_l, parameters)
#     new_conc = append!([0.5], conc_l[2:end])
    conc = log.(conc_l)
    p1 = plot(conc, parameters[1,:], xlabel = "drug conc. [nM]", label="", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "alpha"); ylims!(0.0, 1.2*maximum(parameters[1,:]))

    p2 = plot(conc, parameters[2,:], xlabel = "drug conc. [nM]", label = "", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "beta"); ylims!(0.0, 1.2*maximum(parameters[2,:]))

    p3 = plot(conc, parameters[3,:], xlabel = "drug conc. [nM]", label = "", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "gamma1"); ylims!(0.0, 1.2*maximum(parameters[3,:]))

    p4 = plot(conc, parameters[4,:], xlabel = "drug conc. [nM]", label = "", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "gamma2"); ylims!(0.0, 1.2*maximum(parameters[4,:]))

    plot(p1, p2, p3, p4)
    plot!(size = (800, 400), dpi = 150)
end


# pyplot();
function plotUnitCombin(params::Array, i::Int, title::String, bool::Any, g2::Matrix, g1::Matrix, g2_0::Array, g1_0::Array, concL, ConcGem)
    """ Given estimated parameters for each trial, 
    solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, 
    along with their corresponding real data,
    for a longer time which is 2 times of the 
    original time (~195 hours)
    """
    n_times = range(0.0; stop = 200.0, length = 400)
    alg, n_prob, _ = ddesolve(collect(n_times), g1, g2, g1_0, g2_0, params, i)

    solution = solve(n_prob, alg; constrained=true)

    plot(n_times, solution(n_times, idxs=1).u, label = "G1", xlabel = "time [hours]", ylabel = "number of cells", lw=2.0, alpha = 0.6, color =:green)
    plot!(n_times, solution(n_times, idxs=2).u, label = "G2", legend=bool, legendfontsize=7, fg_legend = :transparent, lw=2.0, alpha = 0.6, color=:sienna)
    plot!(n_times, (solution(n_times, idxs=2).u + solution(n_times, idxs=1).u), label = "total", dpi = 250, lw=2.0, alpha = 0.6, color=:hotpink, margin = 20mm)
    plot!(annotation=(100.5,125.0, text("$concL nM lapat. & $ConcGem nM Dox.", 10)))
end


function plot4combin(ddeparam, g2_l::Matrix, g1_l::Matrix, g2_0_l::Array, g1_0_l::Array, i::Int, conc_l ,conc_g)
    """ here we plot 8 combinations of lapatinib i and all the doxorubicin """ 
    pl = [plotUnitCombin(ddeparam[i,j,:], 1, "", j == 1, g2_l, g1_l, g2_0_l, g1_0_l, conc_l[i], conc_g[j]) for j in 1:8]
    plot(pl..., layout=(2,4))
    plot!(size = (1600, 600), dpi = 250)
    ylims!((0.0, 120.0))
end