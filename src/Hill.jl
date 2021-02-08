"""
This file fits Hill function to the parameters
"""

""" This functions takes in hill parameters for all the concentrations and calculates
DDE parameters, passes them to residual function and based off of these, optimizes the model
and estimates hill parameters. """
function residHill(x::Vector, pControl::Vector, conc::Vector, g1::Matrix, g2::Matrix)

    params = getODEparams(x, conc)
    t = LinRange(0.0, 0.5 * size(g1, 1), size(g1, 1))
    res = 0.0
    # Solve each concentration separately
    for ii = 1:length(conc)
        res += newPredict(params[:, ii], pControl, t, g1[:, ii], g2[:, ii])[1]
    end

    return res
end


""" Generic setup for optimization. """
function optimize_helper(f, low::Vector, high::Vector, maxstep::Int)
    results_ode = bboptimize(
        f;
        SearchRange = collect(zip(low, high)),
        NumDimensions = length(low),
        TraceMode = :verbose,
        TraceInterval = 100,
        MaxSteps = maxstep,
    )

    return best_fitness(results_ode), best_candidate(results_ode)
end


""" Hill optimization function. """
function optimize_hill(conc::Vector, pControl::Vector, g1::Matrix, g2::Matrix; maxstep = 300000)

    f(x) = residHill(x, pControl, conc, g1, g2)

    # [EC50, k, min_a1, max_a1, min_a2, max_a2, min_b1, max_b1, min_b2, max_b2, max_g11, max_g12, max_g21, max_g22, min%G1]
    low = [minimum(conc), 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.25]
    high = [maximum(conc), 10.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.75]

    return optimize_helper(f, low, high, maxstep)
end

""" A function to convert the estimated hill parameters back to ODE parameters. """
function getODEparams(p::Vector, concentrations::Vector{Float64})
    effects = Matrix{eltype(p)}(undef, 9, length(concentrations))

    # Scaled drug effect
    xx = 1.0 ./ (1.0 .+ (p[1] ./ (concentrations .+ eps())) .^ p[2])

    # [EC50, left, right, steepness]
    effects[1, :] = p[3] .+ (p[4] - p[3]) .* xx # a1
    effects[2, :] = p[5] .+ (p[6] - p[5]) .* xx # a2
    effects[3, :] = p[7] .+ (p[8] - p[7]) .* xx # b1
    effects[4, :] = p[9] .+ (p[10] - p[9]) .* xx # b2
    effects[5, :] = p[11] .* xx # death in a1
    effects[6, :] = p[12] .* xx # death in a2
    effects[7, :] = p[13] .* xx # death in b1
    effects[8, :] = p[14] .* xx # death in b2
    effects[9, :] .= p[15] # %G1

    return effects
end
