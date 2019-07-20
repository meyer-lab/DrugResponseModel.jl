module one_state_model

export get_data, optimIt, plotIt, ode_optimIt, ode_plotIt

include("importData.jl")
include("DDEmodel.jl")
include("plot.jl")

end