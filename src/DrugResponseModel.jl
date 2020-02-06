
module DrugResponseModel

using CSV
using DataFrames
using Plots
using Measures
using BlackBoxOptim
using LinearAlgebra
using Base.Threads
using Expokit
import Calculus
import ForwardDiff

include("importData.jl")
include("ODEmodel.jl")
include("Hill.jl")
include("plot.jl")

export setup_data,
    ODEoptimizer,
    ode_plotIt,
    plotIt,
    correlationPlot,
    plot_all,
    ODEplot_all,
    plot_parameters,
    optimize_hill,
    getODEparams,
    sensitivity,
    plotUnitSensitivity,
    allSensitivity,
    ODEplot_allPerc,
    plotGradient,
    predict

end # module
