using Plots, CSV, DataFrames 
include("DDEmodel.jl")
gr()
""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

#-------------------- Plots the long-term predicted behavior given parameters -----------------------#
function plotIt(params::Array, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop::DataFrame, i::Int, title::String, bool::Any)
    """ Given estimated parameters for each trial, 
    solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, 
    along with their corresponding real data,
    for a longer time which is 2 times of the 
    original time (~195 hours)
    """
    lags = [params[3], params[4]]
    t = LinRange(0.0, 95.5, 192)
    t_new = LinRange(0.0, 195.0, 292)
    tspan_new = (0.0, 195.5)
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]
    u0_new = [g1_0[i], g2_0[i]]
    prob_new = DDEProblem(DDEmodel, u0_new, h, tspan_new, params; constant_lags = lags)
    solution = solve(prob_new, MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true))

    plot(t_new, solution(t_new, idxs=1).u, label = "G1 est", dpi = 150, title = title, xlabel = "time [hours]", ylabel = "# of cells", lw=2.0, alpha = 0.6)
    scatter!(t, g1[:, i], label = "G1", dpi = 150, markersize = 1.0, marker=([:dot :d], 1, 0.8, Plots.stroke(0.1, :gray)))
    plot!(t_new, solution(t_new, idxs=2).u, label = "G2 est", legend=bool, legendfontsize=6, fg_legend = :transparent, lw=2.0, alpha = 0.6)
    scatter!(t, g2[:, i], label = "G2", dpi = 150, markersize = 1.0, marker=([:dot :d], 1, 0.8, Plots.stroke(0.1, :gray)))
    plot!(t_new, (solution(t_new, idxs=2).u + solution(t_new, idxs=1).u), label = "total est", dpi = 150, lw=2.0, alpha = 0.6)
    scatter!(t, pop[i], label = "total", dpi = 150, markersize = 1.0, marker=([:dot :d], 1, 0.8, Plots.stroke(0.1, :gray)))
end
