using Plots, CSV, Distributed, DataFrames;
include("DDEmodel.jl")

""" 
        This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot.
"""

#-------------------- Plots the long-term predicted behavior given parameters -----------------------#

function plotIt(params::Array, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop::DataFrame, i::Int, title::String)
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
    h_g1_params, h_g2_params = find_history(g1, g2)
    h(pp, t_new) = [h_g1_params[1]*exp.(t_new*h_g1_params[2]), h_g2_params[1]*exp.(t_new*h_g2_params[2])]
    u0_new = [g1_0[i], g2_0[i]]
    prob_new = DDEProblem(DDEmodel, u0_new, h, tspan_new, params; constant_lags = lags)
    solution = solve(prob_new, MethodOfSteps(Tsit5()))


    plot(t_new, solution(t_new, idxs=1).u, label = "G1 est", dpi = 150, title = title, xlabel = "time [hours]", ylabel = "# of cells")
    plot!(t, g1[:, i], label = "G1", dpi = 150)
    plot!(t_new, solution(t_new, idxs=2).u, label = "G2 est", legend=:topright, dpi = 150)
    plot!(t, g2[:, i], label = "G2", dpi = 150)
    plot!(t_new, (solution(t_new, idxs=2).u + solution(t_new, idxs=1).u), label = "total est", dpi = 150)
    plot!(t, pop[i], label = "total", dpi = 150)
    # savefig("gem_3_dde_long.png")
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

# using CSV

# # Reading the file containing params
# param_lap_dde = CSV.read(".//figures//Lapatinib//params_lap_DDE.csv")
# param_gem_dde = CSV.read(".//figures//Gem//params_gem_DDE.csv")
# param_dox_dde = CSV.read(".//figures//Dox//params_dox_DDE.csv")
# param_taxol1_dde = CSV.read(".//figures//taxol//params_taxol1_DDE.csv")

# # Convert the DataFrame to Matrix for plotting
# lap = convert(Matrix, param_lap_dde[:,2:end])
# gem = convert(Matrix, param_gem_dde[:,2:end])
# dox = convert(Matrix, param_dox_dde[:,2:end])
# tax = convert(Matrix, param_taxol1_dde[:,2:end])
