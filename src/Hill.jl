"""
This file fits Hill function to the parameters
"""

""" This functions takes in hill parameters for all the concentrations and calculates
DDE parameters, passes them to residual function and based off of these, optimizes the model
and estimates hill parameters. """
function residHill(x::Vector, conc::Vector, g1::Matrix, g2::Matrix)

    params = getODEparams(x, conc)
    t = LinRange(0.0, 0.5 * size(g1, 1), size(g1, 1))
    res = 0.0
    # Solve each concentration separately
    for ii = 1:length(conc)
        res += DrugResponseModel.newPredict(params[:, ii, 1], params[:, 1, 1], t, g1[:, ii], g2[:, ii])[1]
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
function optimize_hill(conc::Vector, g1::Matrix, g2::Matrix; maxstep = 300000)

    f(x) = residHill(x, conc, g1, g2)

    # [EC50, k, min_a1, max_a1, min_a2, max_a2, min_b1, max_b1, min_b2, max_b2, max_g11, max_g12, max_g21, max_g22, min%G1]
    low = [minimum(conc), 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.25]
    high = [maximum(conc), 10.0, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0.75]

    return optimize_helper(f, low, high, maxstep)
end

function getODEparams(p, conc)
    if length(p) == 59
        nMax = 5
    elseif length(p) == 15
        nMax = 1
    elseif length(p) == 26
        nMax = 2
    end

    effects = zeros(eltype(p), 9, length(conc[:, 1]), Int((length(p)-4)/11))
    k = 1
    sizep = 11 # the size of independent parameters, meaning except for control.
    j = nMax * sizep + 1 # the starting index of "control parameters", according to the number of drugs being fitted at once.
    # Scaled drug effect
    for i = 1:nMax
        xx = 1.0 ./ (1.0 .+ (p[k] ./ conc[:, i]) .^ p[k + 1])

        # [EC50, left, right, steepness]
        effects[1, :, i] = p[j] .+ (p[k + 2] - p[j]) .* xx
        effects[2, :, i] = p[j + 1] .+ (p[k + 3] - p[j + 1]) .* xx
        effects[3, :, i] = p[j + 2] .+ (p[k + 4] - p[j + 2]) .* xx
        effects[4, :, i] = p[j + 3] .+ (p[k + 5] - p[j + 3]) .* xx
        effects[5, :, i] = p[k + 6] .* xx
        effects[6, :, i] = p[k + 7] .* xx
        effects[7, :, i] = p[k + 8] .* xx
        effects[8, :, i] = p[k + 9] .* xx
        effects[9, :, i] .= p[k + 10]
        
        k += sizep
    end
    return effects
end
