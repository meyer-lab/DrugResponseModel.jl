""" In this file we fit all the drugs att once. """

function getODEparamsAll(p::Array{Float64,1}, concentrations::Array{Float64,2})
    effects = zeros(9, 8, 4)

    k = 1
    # Scaled drug effect
    for i = 1:4
        xx = 1.0 ./ (1.0 .+ (p[k] ./ concentrations[:, i]) .^ p[k + 1])

        effects[1, :, i] = p[25] .+ (p[k + 2] - p[25]) .* xx
        effects[2, :, i] = p[26] .+ (p[k + 3] - p[26]) .* xx
        effects[3, :, i] = p[k + 4] .* xx
        effects[4, :, i] = p[k + 5] .* xx
        k += 6
    end
    effects[5, :, :] .= p[27] #percentage in G1
    effects[6, :, :] .= floor(p[28]) #nG1
    effects[7, :, :] .= floor(p[29]) #nG2
    effects[8, :, :] .= floor(p[30]) #nD1
    effects[9, :, :] .= floor(p[31]) #nD2

    return effects
end

function residHillAll(hillParams::Array{Float64, 1}, concentrations::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3})
    res = Atomic{eltype(hillParams)}(0.0)
    params = getODEparamsAll(hillParams, concentrations)

    # Solve for all drugs
    for j = 1:4
        @threads for ii = 1:length(concentrations[:, j])
            atomic_add!(
                res,
                cost(
                    params[:, ii, j],
                    g1[:, ii, j],
                    g2[:, ii, j],
                    Int(floor(params[6, ii, j])),
                    Int(floor(params[7, ii, j])),
                    Int(floor(params[8, ii, j])),
                    Int(floor(params[9, ii, j])),
                ),
            )
        end
    end

    return res[]
end

""" Hill optimization function for all drugs. """
function optimize_hillAll(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxstep = 1E5)
    hillCostAll(hillParams) = residHillAll(hillParams, concs, g1, g2)

    low = @LArray [minimum(concs[:, 1]), 0.01, 1e-9, 1e-9, 0.0, 0.0, minimum(concs[:, 2]), 0.01, 1e-9, 1e-9, 0.0, 0.0, minimum(concs[:, 3]), 0.01, 1e-9, 1e-9, 0.0, 0.0, minimum(concs[:, 4]), 0.01, 1e-9, 1e-9, 0.0, 0.0, 1e-9, 1e-9, 0.45, 2, 10, 0, 0] (:Lap_EC50, :Lap_steepness, :Lap_maxG1ProgRate, :Lap_maxG2ProgRate, :Lap_maxDeathG1Rate, :Lap_maxDeathG2Rate, :Dox_EC50, :Dox_steepness, :Dox_maxG1ProgRate, :Dox_maxG2ProgRate, :Dox_maxDeathG1Rate, :Dox_maxDeathG2Rate, :Gem_EC50, :Gem_steepness, :Gem_maxG1ProgRate, :Gem_maxG2ProgRate, :Gem_maxDeathG1Rate, :Gem_maxDeathG2Rate, :Tax_EC50, :Tax_steepness, :Tax_maxG1ProgRate, :Tax_maxG2ProgRate, :Tax_maxDeathG1Rate, :Tax_maxDeathG2Rate, :G1ProgRateControl, :G2ProgRateControl, :percG1, :nG1, :nG2, :nD1, :nD2)
    high = [maximum(concs[:, 1]), 10.0, 3.0, 3.0, 1.0, 1.0, maximum(concs[:, 2]), 10.0, 3.0, 3.0, 1.0, 1.0, maximum(concs[:, 3]), 10.0, 3.0, 3.0, 1.0, 1.0, maximum(concs[:, 4]), 10.0, 3.0, 3.0, 1.0, 1.0, 3.0, 3.0, 0.55, 60, 180, 50, 50]

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