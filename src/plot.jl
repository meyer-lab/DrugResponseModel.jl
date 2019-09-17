using Plots, CSV, DataFrames
include("importData.jl")
include("DDEmodel.jl")
pyplot()
""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

#-------------------- Plots the long-term predicted behavior given parameters -----------------------#
function plotIt(params::Array, i, title::String, bool::Any)
    """ Given estimated parameters for each trial, 
    solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, 
    along with their corresponding real data,
    for a longer time which is 2 times of the 
    original time (~195 hours)
    """
    _, pop, g2, g1, g2_0, g1_0 = setup_data("lapatinib")
    lags = [params[3], params[4]]
    times = range(0.0; stop = 95.5, length = 192)
    n_times = range(0.0; stop = 195.5, length = 292)
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]
    u0_new = [g1_0[i], g2_0[i]]
    prob_new = DDEProblem(DDEmodel, u0_new, h, extrema(n_times), params; constant_lags = lags)
    
    solution = solve(prob_new, MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true))

    plot(n_times, solution(n_times, idxs=1).u, label = "G1 est", dpi = 150, xlabel = "time [hours]", ylabel = "# of cells", lw=2.0, alpha = 0.6, color =:green)
    plot!(times, g1[:, i], label = "G1", dpi = 150, markersize = 1.0, color=:darkgreen)
    plot!(n_times, solution(n_times, idxs=2).u, label = "G2 est", legend=bool, legendfontsize=7, fg_legend = :transparent, lw=2.0, alpha = 0.6, color=:sienna)
    plot!(times, g2[:, i], label = "G2", dpi = 150, markersize = 1.0, color=:darkorange)
    plot!(n_times, (solution(n_times, idxs=2).u + solution(n_times, idxs=1).u), label = "total est", dpi = 150, lw=2.0, alpha = 0.6, color=:hotpink)
    plot!(times, pop[i], label = "total", dpi = 150, markersize = 1.0, color=:indigo)
    plot!( annotation=[ (75,90, text(title, 12)) ])
end
