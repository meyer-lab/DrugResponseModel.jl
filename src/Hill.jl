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
    )

    return best_fitness(results_ode), best_candidate(results_ode)
end


""" Hill optimization function. """
function optimize_hill(conc::Vector, g1::Matrix, g2::Matrix; maxstep = 300000)

    f(x) = residHill(x, conc, g1, g2)

    # [EC50, k, max_a1, max_a2, max_b1,  max_b2, max_b3, max_b4, max_g11, max_g12, max_g21, max_g22, max_g23, max_g24, min_a1, min_a2, min_b1,  min_b2, min_b3, min_b4]
    low = [minimum(conc); 1e-9 * ones(25)]
    high = [maximum(conc); 10.0; 2.5 * ones(24)]

    return optimize_helper(f, low, high, maxstep)
end

# function optimC(g1, g2)

#     ff(x) = residC(x, g1, g2)
#     low = append!(1e-6 * ones(6), zeros(6))
#     high = append!(3.0 * ones(6), zeros(6))
#     return optimize_helper(ff, low, high, 300000)
# end

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
        effects[10, :, i] = p[k + 11] .* xx  # g12
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

function find_cellnumber_ec50()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    _, popul2, g1s2, g2s2 = load(189, 2)
    _, popul3, g1s3, g2s3 = load(189, 3)

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims = 4)
    g2S = cat(g2s1, g2s2, g2s3, dims = 4)
    g1m = mean(g1S, dims = 4) # mean G1
    g2m = mean(g2S, dims = 4) # mean G2

    total = g1m[189, :, :, 1] + g2m[189, :, :, 1]

    f(x) = costingss(x, total, concs)

    low = [5.0, 0.01, 1e-9, 1e-9, 5.0, 0.01, 1e-9, 1e-9, 1.0, 0.01, 1e-9, 1e-9, 0.5, 0.01, 1e-9, 1e-9, 5.0, 0.01, 1e-9, 1e-9]
    high = [500.0, 10.0, 3.0, 3.0, 500.0, 10.0, 3.0, 3.0, 200.0, 10.0, 3.0, 3.0, 100.0, 10.0, 3.0, 3.0, 500.0, 10.0, 3.0, 3.0]
    _, p = optimize_helper(f, low, high, 100000)

    k = 1
    num = zeros(8, 5)
    for i = 1:5
        num[:, i] = hill_func(p[k:(k + 3)], concs[:, i])
        k += 4
    end
    p1 = plot(concs[:, 1], num[:, 1], label = "model", title = "lapatinib")
    plot!(concs[:, 1], total[:, 1], label = "data")
    p2 = plot(concs[:, 2], num[:, 2], label = "model", title = "doxorubicin")
    plot!(concs[:, 2], total[:, 2], label = "data")
    p3 = plot(concs[:, 3], num[:, 3], label = "model", title = "gemcitabine")
    plot!(concs[:, 3], total[:, 3], label = "data")
    p4 = plot(concs[:, 4], num[:, 4], label = "model", title = "paclitaxel")
    plot!(concs[:, 4], total[:, 4], label = "data")
    p5 = plot(concs[:, 5], num[:, 5], label = "model", title = "palbociclib")
    plot!(concs[:, 5], total[:, 5], label = "data")
    plt = plot(p1, p2, p3, p4, p5)
    savefig(plt, "hils.svg")
    optimized_params = [
        57.3004,
        2.00363,
        0.948318,
        3.0,
        11.2302,
        2.01893,
        0.258156,
        3.0,
        8.95294,
        3.16891,
        0.395267,
        3.0,
        2.46285,
        4.53685,
        0.274721,
        3.0,
        30.1616,
        1.95144,
        1.38786,
        3.0,
    ]
end
