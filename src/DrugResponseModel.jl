
module DrugResponseModel

using CSV
using DataFrames
using Plots
using Measures
using BlackBoxOptim
using LinearAlgebra
using Base.Threads
import Calculus

include("importData.jl")
include("ODEmodel.jl")
include("Hill.jl")
include("plot.jl")
include("allDrugs.jl")

export setup_data,
    load,
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
    plotGradient

end # module
