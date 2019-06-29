""" This file contains a function to plot the parameters against the  drug concentrations, in a 2x2 subplot."""
using Plots, CSV;

# Reading the file containing params
param_lap_dde = CSV.read(".//figures//Lapatinib//dde//params_lap_DDE.csv")
param_gem_dde = CSV.read(".//figures//Gem//dde//params_gem_DDE.csv")
param_dox_dde = CSV.read(".//figures//Dox//dde//params_dox_DDE.csv")
param_taxol1_dde = CSV.read(".//figures//taxol//dde//params_taxol1_DDE.csv")

# Convert the DataFrame to Matrix for plotting
lap = convert(Matrix, param_lap_dde[:,3:end])
gem = convert(Matrix, param_gem_dde[:,3:end])
dox = convert(Matrix, param_dox_dde[:,3:end])
tax = convert(Matrix, param_taxol1_dde[:,3:end])


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
    p1 = plot(gem[8, :], gem[i, :], label = "Gemcitabine", title = param, xlabel = "drug conc. [nM]", ylabel = "param", yticks = 0.0:0.02:0.5)
    p2 = plot(lap[8, :], lap[i, :], label = "Lapatinib", xlabel = "drug conc. [nM]", ylabel = "param", yticks = 0.0:0.05:0.5)
    p3 = plot(dox[8, :], dox[i, :], label = "Doxorubicin", xlabel = "drug conc. [nM]", ylabel = "param", yticks = 0.0:0.05:0.5)
    p4 = plot(tax[8, :], tax[i, :], label = "Paclitaxel", xlabel = "drug conc. [nM]", ylabel = "param", yticks = 0.0:0.1:1.0)
    plot(p1, p2, p3, p4, dpi = 100)
end

plot_param_conc(lap, gem, dox, tax, 1, "alpha")
# We can save the figure by using the following line with an appropriate name in the paramthesis.
# savefig("")