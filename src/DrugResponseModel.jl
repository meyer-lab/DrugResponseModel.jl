module DrugResponseModel

export get_data, remove_peaks, ddesolve, optimization, PlotIt, residHill, optimize_hill

include("DDEmodel.jl")
include("ODEmodel.jl")
include("plot.jl")
include("Hill.jl")

end # module
