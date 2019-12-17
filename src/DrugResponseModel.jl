
__precompile__()
module DrugResponseModel

using CSV
using DataFrames
using Plots
using Measures
using BlackBoxOptim
using OrdinaryDiffEq

include("importData.jl")
include("ODE4.jl")
include("ODEmodel.jl")
include("plot.jl")

export setup_data, remove_peaks, ODEoptimizer, ode_plotIt, plotIt, correlationPlot, plot_all, plot_parameters, ODEplot_all, ODEoptimizer4, ode_plotIt4, update_coef4, ODEmodel4, cost, predict, ODEplot_all4

end # module
