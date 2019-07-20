module one_state_model

export get_data, DDEmodel, DDEsolve, resid, optimIt, ODEmodel, ODEsolve, residuals, ode_optimIt, ode_plotIt

include("importData.jl")
include("DDEmodel.jl")
include("plot.jl")

end