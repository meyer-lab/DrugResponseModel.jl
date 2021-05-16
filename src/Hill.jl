"""
This file fits Hill function to the parameters
"""

""" This functions takes in hill parameters for all the concentrations and calculates
DDE parameters, passes them to residual function and based off of these, optimizes the model
and estimates hill parameters. """
function residHill(x::Vector, conc::Vector, g1::Matrix, g2::Matrix)

    res = 0.0
    for i=3:8
        res += 60*(maximum([0, (x[i] - x[i + 12])]))^2
    end
    params = getODEparams(x, conc)
    t = LinRange(0.0, 0.5 * size(g1, 1), size(g1, 1))
    # Solve each concentration separately
    for ii = 1:length(conc)
        res += predict(params[:, ii, 1], params[:, 1, 1], t, g1[:, ii], g2[:, ii])[1]
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

    # [EC50, k, max_a1, max_a2, max_b1,  max_b2, max_b3, max_b4, max_g11, max_g12, max_g21, max_g22, max_g23, max_g24, min_a1, min_a2, min_b1,  min_b2, min_b3, min_b4]
    low = [minimum(conc); 1e-9 * ones(19)]
    high = [maximum(conc); 10.0; 2.5 * ones(18)]

    return optimize_helper(f, low, high, maxstep)
end

function getODEparams(p, concs)
    if size(concs) == (8,)
        conc = zeros(8, 1)
        conc[:, 1] = concs
    else
        conc = concs
    end

    nMax = Int((length(p) - 6) / 14)

    effects = zeros(eltype(p), 12, length(conc[:, 1]), nMax)
    k = 1
    sizep = 14 # the size of independent parameters, meaning except for control.
    j = nMax * sizep + 1 # the starting index of "control parameters", according to the number of drugs being fitted at once.
    # Scaled drug effect
    for i = 1:nMax
        xx = 1.0 ./ (1.0 .+ (p[k] ./ conc[:, i]) .^ p[k + 1])

        # [EC50, left, right, steepness]
        effects[1, :, i] = p[j] .+ (p[k + 2] - p[j]) .* xx # a1
        effects[2, :, i] = p[j + 1] .+ (p[k + 3] - p[j + 1]) .* xx # a2
        effects[3, :, i] = p[j + 2] .+ (p[k + 4] - p[j + 2]) .* xx # b1
        effects[4, :, i] = p[j + 3] .+ (p[k + 5] - p[j + 3]) .* xx # b2
        effects[5, :, i] = p[j + 4] .+ (p[k + 6] - p[j + 4]) .* xx # b3
        effects[6, :, i] = p[j + 5] .+ (p[k + 7] - p[j + 5]) .* xx # b4
        effects[7, :, i] = p[k + 8] .* xx   # g11
        effects[8, :, i] = p[k + 9] .* xx   # g12
        effects[9, :, i] = p[k + 10] .* xx  # g21
        effects[10, :, i] = p[k + 11] .* xx # g22
        effects[11, :, i] = p[k + 12] .* xx # g23
        effects[12, :, i] = p[k + 13] .* xx # g24

        k += sizep
    end
    return effects
end
