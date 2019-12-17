
__precompile__()
module DrugResponseModel

using CSV
using DataFrames
using Plots
using Measures
using BlackBoxOptim
import ForwardDiff
import LinearAlgebra

include("importData.jl")
include("ODEmodel.jl")
include("plot.jl")

export setup_data, ODEoptimizer, ode_plotIt, plotIt, correlationPlot, plot_all, ODEplot_all

end # module
