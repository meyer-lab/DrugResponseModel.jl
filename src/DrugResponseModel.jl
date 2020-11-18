
module DrugResponseModel

ENV["GKSwstype"] = "100"

using CSV
using DataFrames
using Plots
using Measures
using LinearAlgebra
using Statistics
using DSP: conv
using Optim
import ForwardDiff
import ExponentialUtilities
import LineSearches
using Polynomials

include("importData.jl")
include("ODEmodel.jl")
include("Hill.jl")
include("plot.jl")
include("allDrugs.jl")
include("combination.jl")
include("replicates.jl")

export setup_data,
    load,
    ode_plotIt,
    plotIt,
    correlationPlot,
    plot_all,
    optimize_hill,
    getODEparams,
    ODEplot_allPerc,
    numcells,
    getODEparamsAll,
    numcells,
    CombinationParam,
    fullCombinationParam,
    plotEffectsCombin,
    plotNumcells,
    blissCellNum,
    predict,
    savitzky_golay_filter,
    optim_all,
    BlissModelComb,
    Heatmap

end # module
