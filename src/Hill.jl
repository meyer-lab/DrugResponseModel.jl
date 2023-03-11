"""
This file fits Hill function to the parameters
"""

""" This functions takes in hill parameters for all the concentrations and calculates
DDE parameters, passes them to residual function and based off of these, optimizes the model
and estimates hill parameters. """
function residHill(x::Vector, conc::Vector, g1::Matrix, g2::Matrix)

    res = 0.0
    for i = 3:10
        res += 40 * (maximum([0, (x[i] - x[i + 16])]))^2
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
        Method = :adaptive_de_rand_1_bin_radiuslimited
    )

    return best_fitness(results_ode), best_candidate(results_ode)
end


""" Hill optimization function. """
function optimize_hill(conc::Vector, g1::Matrix, g2::Matrix; maxstep = 200000)

    f(x) = residHill(x, conc, g1, g2)

    # [EC50, k, max_a1, max_a2, max_b1,  max_b2, max_b3, max_b4, max_g11, max_g12, max_g21, max_g22, max_g23, max_g24, min_a1, min_a2, min_b1,  min_b2, min_b3, min_b4]
    low = [minimum(conc); 1e-9 * ones(25)]
    high = [2*maximum(conc); 50.0; 4.0 * ones(24)]

    return optimize_helper(f, low, high, maxstep)
end


function getODEparams(p, conc)
    nMax = Int((length(p) - 8) / 18)

    effects = zeros(eltype(p), 16, length(conc[:, 1]), nMax)
    k = 1
    sizep = 18 # the size of independent parameters, meaning except for control.
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
        effects[7, :, i] = p[j + 6] .+ (p[k + 8] - p[j + 6]) .* xx # b3
        effects[8, :, i] = p[j + 7] .+ (p[k + 9] - p[j + 7]) .* xx # b4
        effects[9, :, i] = p[k + 10] .* xx   # g11
        effects[10, :, i] = p[k + 11] .* xx   # g12
        effects[11, :, i] = p[k + 12] .* xx  # g21
        effects[12, :, i] = p[k + 13] .* xx # g22
        effects[13, :, i] = p[k + 14] .* xx # g23
        effects[14, :, i] = p[k + 15] .* xx # g24
        effects[15, :, i] = p[k + 16] .* xx # g23
        effects[16, :, i] = p[k + 17] .* xx # g24

        k += sizep
    end
    return effects
end

hill_func(p, conc) = p[4] .+ (p[3] - p[4]) ./ (1.0 .+ (p[1] ./ conc) .^ p[2])


function costingss(pp, total, concs)
    cost = 0
    k = 1
    for i = 1:5
        cost += sum((hill_func(pp[k:(k + 3)], concs[:, i]) .- total[:, i]) .^ 2)
        k += 4
    end
    return cost
end
