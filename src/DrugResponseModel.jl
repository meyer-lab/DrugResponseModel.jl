
module DrugResponseModel

ENV["GKSwstype"] = "100"

using CSV
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
include("combination.jl")

export setup_data,
    load,
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
    numcells,
    getODEparamsAll,
    optimize_hillAll,
    BlissCombination,
    fullCombinationParam,
    plotEffectsCombin,
    plotNumcells,
    combin2drugs,
    blissCellNum,
    temporal_combination,
    plotTemporalCombin,
    find_IC50

end # module
