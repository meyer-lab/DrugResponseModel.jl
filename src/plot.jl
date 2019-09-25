using Plots, CSV, DataFrames
include("importData.jl")
include("DDEmodel.jl")
gr()
""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

#-------------------- Plots the long-term predicted behavior given parameters -----------------------#
function plotIt(params::Array, i::Int, title::String, bool::Any)
    """ Given estimated parameters for each trial, 
    solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, 
    along with their corresponding real data,
    for a longer time which is 2 times of the 
    original time (~195 hours)
    """
    _, pop, g2, g1, g2_0, g1_0 = setup_data("lapatinib")
    times = range(0.0; stop = 95.5, length = 192)
    n_times = range(0.0; stop = 200.0, length = 400)
    alg, n_prob, _ = ddesolve(n_times, g1, g2, g1_0, g2_0, params, i)

    solution = solve(n_prob, alg; constrained=true)

    plot(n_times, solution(n_times, idxs=1).u, label = "G1 est", dpi = 150, xlabel = "time [hours]", ylabel = "# of cells", lw=2.0, alpha = 0.6, color =:green)
    plot!(times, g1[:, i], label = "G1", dpi = 150, markersize = 1.0, color=:darkgreen)
    plot!(n_times, solution(n_times, idxs=2).u, label = "G2 est", legend=bool, legendfontsize=7, fg_legend = :transparent, lw=2.0, alpha = 0.6, color=:sienna)
    plot!(times, g2[:, i], label = "G2", dpi = 150, markersize = 1.0, color=:darkorange)
    plot!(n_times, (solution(n_times, idxs=2).u + solution(n_times, idxs=1).u), label = "total est", dpi = 150, lw=2.0, alpha = 0.6, color=:hotpink)
    plot!(times, pop[i], label = "total", dpi = 150, markersize = 1.0, color=:indigo)
    plot!( annotation=[ (75,90, text(title, 12)) ])
end

#-------------------------- plot for the G1 G2 correlation ---------------------------#
function correlationPlot(g1, g2, labels, xlabel, ylabel, ymax)

    p1 = scatter(g1[:,1], g2[:,1], title = labels[1], xlabel = xlabel, ylabel=ylabel, alpha = 0.6, color=[:gray], line=:dot, marker=(:dot, 1.5))
    p2 = scatter(g1[:,2], g2[:,2], title = labels[2], xlabel = xlabel, ylabel=ylabel, alpha = 0.6, color=[:gray], line=:dot, marker=(:dot, 1.5))
    p3 = scatter(g1[:,3], g2[:,3], title = labels[3], xlabel = xlabel, ylabel=ylabel, alpha = 0.6, color=[:gray], line=:dot, marker=(:dot, 1.5))
    p4 = scatter(g1[:,4], g2[:,4], title = labels[4], xlabel = xlabel, ylabel=ylabel, alpha = 0.6, color=[:gray], line=:dot, marker=(:dot, 1.5))
    p5 = scatter(g1[:,5], g2[:,5], title = labels[5], xlabel = xlabel, ylabel=ylabel, alpha = 0.6, color=[:gray], line=:dot, marker=(:dot, 1.5))
    p6 = scatter(g1[:,6], g2[:,6], title = labels[6], xlabel = xlabel, ylabel=ylabel, alpha = 0.6, color=[:gray], line=:dot, marker=(:dot, 1.5))
    p7 = scatter(g1[:,7], g2[:,7], title = labels[7], xlabel = xlabel, ylabel=ylabel, alpha = 0.6, color=[:gray], line=:dot, marker=(:dot, 1.5))
    p8 = scatter(g1[:,8], g2[:,8], title = labels[8], xlabel = xlabel, ylabel=ylabel, alpha = 0.6, color=[:gray], line=:dot, marker=(:dot, 1.5))
    plot(p1, p2, p3, p4, p5, p6, p7, p8, legend=:false, layout=(2,4))
    plot!(size=(1200,600), dpi=200)
    ylims!((0, ymax))
    xlims!((0, ymax))
    savefig("corr.png")
end