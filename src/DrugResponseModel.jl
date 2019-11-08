
__precompile__()
module DrugResponseModel

using DelayDiffEq, DiffEqParamEstim, DataFrames, BlackBoxOptim, Plots, CSV, DataFrames

include("ODEmodel.jl")
include("plot.jl")
include("Hill.jl")


export setup_data, remove_peaks, find_history, ODEoptimizer, ode_plotIt, ddesolve, optimization, plotIt, residHill, optimize_hill, getDDEparams, ParamForBliss, BlissCombination, correlationPlot, plot_all, plot_parameters, ODEplot_all, plot3D

precompile(optimization, (Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Int, Vector{Float64}, Vector{Float64}, Int))

end # module
