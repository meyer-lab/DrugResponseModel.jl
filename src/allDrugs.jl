""" In this file we fit all the drugs att once. """

function getODEparamsAll(p::Array{Float64, 1}, concentrations::Array{Float64, 2})
    effects = zeros(9, length(concentrations[:, 1]), 5)

    k = 1
    # Scaled drug effect
    for i = 1:5
        xx = 1.0 ./ (1.0 .+ (p[k] ./ concentrations[:, i]) .^ p[k + 1])

        effects[1, :, i] = p[36] .+ (p[k + 2] - p[36]) .* xx
        effects[2, :, i] = p[37] .+ (p[k + 3] - p[37]) .* xx
        effects[3, :, i] = p[k + 4] .* xx
        effects[4, :, i] = p[k + 5] .* xx
        effects[5, :, i] .= p[k + 6] #percentage in G1
        k += 7
    end

    effects[6, :, :] .= 1
    effects[7, :, :] .= 1
    effects[8, :, :] .= 0
    effects[9, :, :] .= 0

    return effects
end

function residHillAll(hillParams::Vector, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0

    # Solve for all drugs
    t = 1
    for j = 1:5
        hill = [
            hillParams[t],
            hillParams[36],
            hillParams[t + 2],
            hillParams[t + 1],
            hillParams[37],
            hillParams[t + 3],
            hillParams[t + 4],
            hillParams[t + 5],
            hillParams[t + 6],
            1,
            1,
            0,
            0,
        ]
        t += 7
        res += residHill(hill, concentrations[:, j], g1[:, :, j], g2[:, :, j])
    end

    return res
end

""" Hill optimization function for all drugs. """
function optimize_hillAll(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxstep = 1E5)
    hillCostAll(hillParams) = residHillAll(hillParams, concs, g1, g2)

    # The parameters used here in order:
    #(:Lap_EC50, :Lap_steepness, :Lap_maxG1ProgRate, :Lap_maxG2ProgRate, :Lap_maxDeathG1Rate, :Lap_maxDeathG2Rate, :Dox_EC50, :Dox_steepness, :Dox_maxG1ProgRate, :Dox_maxG2ProgRate, :Dox_maxDeathG1Rate, :Dox_maxDeathG2Rate, :Gem_EC50, :Gem_steepness, :Gem_maxG1ProgRate, :Gem_maxG2ProgRate, :Gem_maxDeathG1Rate, :Gem_maxDeathG2Rate, :Tax_EC50, :Tax_steepness, :Tax_maxG1ProgRate, :Tax_maxG2ProgRate, :Tax_maxDeathG1Rate, :Tax_maxDeathG2Rate, :pal_EC50, :pal_steepness, :pal_maxG1ProgRate, :pal_maxG2ProgRate, :pal_maxDeathG1Rate, :pal_maxDeathG2Rate, :G1ProgRateControl, :G2ProgRateControl, :percG1, :nG1, :nG2, :nD1, :nD2)
    lowPiece = [0.01, 0.05, 0.05, 0.0, 0.0, 0.3]
    low = vcat(
        minimum(concs[:, 1]),
        lowPiece,
        minimum(concs[:, 2]),
        lowPiece,
        minimum(concs[:, 3]),
        lowPiece,
        minimum(concs[:, 4]),
        lowPiece,
        minimum(concs[:, 5]),
        lowPiece,
        1e-9,
        1e-9,
        3,
        3,
        0,
        0,
    )
    highPiece = [10.0, 5.0, 5.0, 1.0, 1.0, 0.7]
    high = vcat(
        maximum(concs[:, 1]),
        highPiece,
        maximum(concs[:, 2]),
        highPiece,
        maximum(concs[:, 3]),
        highPiece,
        maximum(concs[:, 4]),
        highPiece,
        maximum(concs[:, 5]),
        highPiece,
        3.0,
        3.0,
        50,
        70,
        50,
        50,
    )

    results_ode = bboptimize(
        hillCostAll;
        SearchRange = collect(zip(low, high)),
        NumDimensions = length(low),
        TraceMode = :verbose,
        TraceInterval = 100,
        MaxSteps = maxstep,
    )

    return best_fitness(results_ode), best_candidate(results_ode)
end

function optim_all(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}, initial_x)
    f(hillParams) = residHillAll(hillParams, concs, g1, g2)

    function g!(G, hillParams)
        ForwardDiff.gradient!(G, f, hillParams)
        costCenter = f(hillParams)

        # Handle the integer-valued parameters
        for ii = 34:37
            pp = copy(hillParams)
            pp[ii] += 1.0
            G[ii] = f(pp) - costCenter
        end
    end

    low = vcat(ones(32) * 1.0e-9, 0.4, 3, 10, 0, 0)
    hP = [1000.0, 10.0, 7.0, 3.0, 3.0, 3.0]
    high = vcat(hP, hP, hP, hP, hP, 3.0, 3.0, 0.6, 10, 25, 50, 50)

    initial_x = [
        488.5,
        1.2,
        3.0,
        2.3,
        1.1,
        0.47,
        499.8,
        1.1,
        3.0,
        1.6,
        2.1,
        2.4,
        99.9,
        1.3,
        3.0,
        2.3,
        0.42,
        0.97,
        3.8,
        3.3,
        3.0,
        1.7,
        0.89,
        0.40,
        499.5,
        0.86,
        3.0,
        2.3,
        0.94,
        0.35,
        0.37,
        0.44,
        0.55,
        9.4,
        15.4,
        49.8,
        12.5,
    ]

    options = Optim.Options(outer_iterations = 2, show_trace = true, iterations = 3)
    results = optimize(f, g!, low, high, initial_x, Fminbox(GradientDescent()), options)

    return Optim.minimizer(results)
end

# for plotting all the concentrations on top of each other, from average of three reps.
# function avgAllConcs(g0, g1m, g2m, p)
    