
__precompile__()
module DrugResponseModel

using CSV
using DataFrames
using Plots
using Measures
using LsqFit
using DelayDiffEq
using DiffEqParamEstim
using BlackBoxOptim
using OrdinaryDiffEq
using DiffEqOperators
using DiffEqBase

include("importData.jl")
include("ODEmodel.jl")
include("ODE4.jl")
include("plot.jl")
include("DDEmodel.jl")
include("Hill.jl")

export setup_data, remove_peaks, find_history, ODEoptimizer, ode_plotIt, ddesolve, optimization, plotIt, residHill, optimize_hill, getDDEparams, ParamForBliss, BlissCombination, correlationPlot, plot_all, plot_parameters, ODEplot_all, ODEoptimizer4, ode_plotIt4, update_coef4, ODEmodel4, cost, predict, ODEplot_all4

precompile(optimization, (Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Int, Vector{Float64}, Vector{Float64}, Int))

end # module
