"""
This file fits Hill function to the parameters.
"""

""" This functions takes in hill parameters for all the concentrations and calculates
DDE parameters, passes them to residual function and based off of these, optimizes the model
and estimates hill parameters. """
function residHill(
    hillParams::Vector,
    concentrations::Vector{Float64},
    g1::Matrix{Float64},
    g2::Matrix{Float64},
    g1_0::Vector{Float64},
    g2_0::Vector{Float64},
)::Float64
    res = Atomic{eltype(hillParams)}(0.0)
    params = getODEparams(hillParams, concentrations)

    # Solve each concentration separately
    @threads for ii = 1:length(concentrations)
        atomic_add!(res, cost(params[1:5, ii], g1_0[ii], g2_0[ii], g1[:, ii], g2[:, ii], Int(floor(params[6, ii])), Int(floor(params[7, ii]))))
    end

    return res[]
end

""" Hill optimization function. """
function optimize_hill(
    conc_l::Array{Float64, 1},
    g1::Array{Float64, 2},
    g2::Array{Float64, 2},
    g1_0::Array{Float64, 1},
    g2_0::Array{Float64, 1};
    maxstep = 1E5,
)
    hillCost(hillParams) = residHill(hillParams, conc_l, g1, g2, g1_0, g2_0)

    low = [minimum(conc_l), 1e-9, 1e-9, 0.1, 1e-9, 1e-9, 0.0, 0.0, 0.45, 2, 10]
    high = [maximum(conc_l), 3.0, 3.0, 10.0, 3.0, 3.0, 1.0, 1.0, 0.55, 60, 180]

    results_ode = bboptimize(
        hillCost;
        SearchRange = collect(zip(low, high)),
        NumDimensions = length(low),
        TraceMode = :verbose,
        TraceInterval = 100,
        MaxSteps = maxstep,
    )

    return best_fitness(results_ode), best_candidate(results_ode)
end

""" A function to convert the estimated hill parameters back to ODE parameters. """
function getODEparams(p::Vector, concentrations::Vector{Float64})
    effects = Matrix{eltype(p)}(undef, 7, 8)

    # Scaled drug effect
    xx = 1.0 ./ (1.0 .+ (p[1] ./ concentrations) .^ p[4])

    # [EC50, left, right, steepness]
    effects[1, :] = p[2] .+ (p[3] - p[2]) .* xx
    effects[2, :] = p[5] .+ (p[6] - p[5]) .* xx
    effects[3, :] = p[7] .* xx
    effects[4, :] = p[8] .* xx
    effects[5, :] .= p[9]
    effects[6, :] .= floor(p[10])
    effects[7, :] .= floor(p[11])

    return effects
end

""" To find the sensitivity of the model to parameters. """
function sensitivity(
    params::Vector,
    conc::Vector{Float64},
    i::Int,
    g1_0::Vector{Float64},
    g2_0::Vector{Float64},
    g1::Matrix{Float64},
    g2::Matrix{Float64}
)
    result = zeros(200)
    paramRange = LinRange(0.1*params[i], 10*params[i], 200)
    for j=1:200
        params[i]= paramRange[j]
        result[j] = residHill(params, conc, g1, g2, g1_0, g2_0)
    end
    return result, paramRange
end
