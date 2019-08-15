using Plots, CSV, Distributed, DataFrames;
include("DDEmodel.jl")

""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

#-------------------- Plots the long-term predicted behavior given parameters -----------------------#

function PlotIt(g1_0, g2_0, j, min_p, data)
    times = range(0.0; stop = 95.5, length = 192)
    new_times = range(0.0; stop=200.0, length=400)
    tspan_new = (0.0, 200.0)
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()))
    prob_new = DDEProblem(DDEmodel, [g1_0[j], g2_0[j]], h, tspan_new, min_p; constant_lags = [min_p[3], min_p[4]])
    solution = solve(prob_new, alg)
    plot(times, data', show = true, label = ["G1", "G2"], xlabel="time[hours]", ylabel="# of cells", lw=2, legend=:best)
    plot!(new_times, solution(new_times, idxs=1).u, label = "G1 estimated",lw=2)
    plot!(new_times, solution(new_times, idxs=2).u, label = "G2 estimated", lw=2)
end

#--------------------- Plot parameteres versus drug concentration ------------------------#

function plot_param_conc(lap, gem, dox, tax, i, param)
    """ This function to plot parameter vs. concentraition.
    Arguments:
    ----------
    lap, gem, dox, tax: DataFrame containing parameters of treating AU565 human breast cancer cells with Lapatinib, Gemcitabine, Doxorubicin, and   Paclitaxel1, respectively.
    i {int.}: If it is for DDE: A number between 1 and 7, which points to the parameters: [alpha, beta, tau1, tau2, history, gamma1, gamma2]
              If it is for ODE: A number between 1 and 4, which points to the parameters: [alpha, beta, gamma1, gamma2]
    param {str.}: a string used for the title of plots, to show which parameter we are plotting.
    Returns:
    --------
    Returns a 2x2 plot for four drugs. 
    
    """
    p1 = plot(gem[8, :], gem[i, :], label = "Gemcitabine", title = param, xlabel = "drug conc. [nM]", ylabel = "param", yticks = 0.0:round(maximum(gem[i, :])/5 ,digits = 3):maximum(gem[i, :]))
    p2 = plot(lap[8, :], lap[i, :], label = "Lapatinib", xlabel = "drug conc. [nM]", ylabel = "param", yticks = 0.0:round(maximum(lap[i, :])/5, digits = 3):maximum(lap[i, :]))
    p3 = plot(dox[8, :], dox[i, :], label = "Doxorubicin", xlabel = "drug conc. [nM]", ylabel = "param", yticks = 0.0:round(maximum(dox[i, :])/5, digits = 3):maximum(dox[i, :]))
    p4 = plot(tax[8, :], tax[i, :], label = "Paclitaxel", xlabel = "drug conc. [nM]", ylabel = "param", yticks = 0.0:round(maximum(tax[i, :])/5, digits = 3):maximum(tax[i, :]))
    plot(p1, p2, p3, p4, dpi = 100)
end
