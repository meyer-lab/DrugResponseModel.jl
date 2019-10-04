using Plots, CSV, DataFrames
include("importData.jl")
include("DDEmodel.jl")
gr()
""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

#-------------------- Plots the long-term predicted behavior given parameters -----------------------#
function plotIt(params::Array, i::Int, title::String, bool::Any, pop, g2::Matrix, g1::Matrix, g2_0::Array, g1_0::Array)
    """ Given estimated parameters for each trial, 
    solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, 
    along with their corresponding real data,
    for a longer time which is 2 times of the 
    original time (~195 hours)
    """
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
function correlationPlot(g1::Matrix, g2::Matrix, labels::Array, xlabel::String, ylabel::String, ymax::Int)

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
#     savefig("corr.png")
end


function plot_all(parameters, pop, g2::Matrix, g1::Matrix, g2_0::Array, g1_0::Array)
    # i showas the trial number, which could be from 1:control, ..., 8: maximum drug concentraation
    p1 = plotIt(parameters[:, 1], 1, "", false, pop, g2, g1, g2_0, g1_0)
    p2 = plotIt(parameters[:, 2], 2, "", false, pop, g2, g1, g2_0, g1_0)
    p3 = plotIt(parameters[:, 3], 3, "", false, pop, g2, g1, g2_0, g1_0)
    p4 = plotIt(parameters[:, 4], 4, "", false, pop, g2, g1, g2_0, g1_0)
    p5 = plotIt(parameters[:, 5], 5, "", false, pop, g2, g1, g2_0, g1_0)
    p6 = plotIt(parameters[:, 6], 6, "", false, pop, g2, g1, g2_0, g1_0)
    p7 = plotIt(parameters[:, 7], 7, "", false, pop, g2, g1, g2_0, g1_0)
    p8 = plotIt(parameters[:, 8], 8, "", :topleft, pop, g2, g1, g2_0, g1_0)
    plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(2,4))
    plot!(size = (1200, 600), dpi = 200)
    ylims!((0.0, 120.0))
end

function plot_parameters(conc_l::Array, parameters)
    new_conc = append!([0.5], conc_l[2:end])
    conc = log.(new_conc)
    p1 = plot(conc, parameters[1,:], xlabel = "drug conc. [nM]", label="", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "alpha"); ylims!(0.0, 1.2*maximum(parameters[1,:]))

    p2 = plot(conc, parameters[2,:], xlabel = "drug conc. [nM]", label = "", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "beta"); ylims!(0.0, 1.2*maximum(parameters[2,:]))

    p3 = plot(conc, parameters[3,:], xlabel = "drug conc. [nM]", label="", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "tau1"); ylims!(0.0, 1.2*maximum(parameters[3,:]))

    p4 = plot(conc, parameters[4,:], xlabel = "drug conc. [nM]", label = "", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "tau2"); ylims!(0.0, 1.2*maximum(parameters[4,:]))

    p5 = plot(conc, parameters[5,:], xlabel = "drug conc. [nM]", label = "", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "gamma1"); ylims!(0.0, 1.2*maximum(parameters[5,:]))

    p6 = plot(conc, parameters[6,:], xlabel = "drug conc. [nM]", label = "", lw= 2.0, alpha = 0.6, color=[:black :gray], line=(:dot, 1), marker=([:dot :d], 3, 0.7, stroke(0.1, 0.6, :gray)),
        ylabel = "gamma2"); ylims!(0.0, 1.2*maximum(parameters[6,:]))

    plot(p1, p2, p3, p4, p5, p6)
    plot!(size = (1200, 600), dpi = 200)
end
