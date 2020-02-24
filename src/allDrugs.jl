""" In this file we fit all the drugs att once. """

function getODEparamsAll(p::Array{Float64,1}, concentrations::Array{Float64,2})
    effects = effects = zeros(9, 8, 4)

    k = 1
    # Scaled drug effect
    for i=1:4
        xx = 1.0 ./ (1.0 .+ (p[k] ./ concentrations[:, i]) .^ p[k+1])

        effects[1, :, i] = p[k+2] .+ (p[k+3] - p[k+2]) .* xx
        effects[2, :, i] = p[k+4] .+ (p[k+5] - p[k+4]) .* xx
        effects[3, :, i] = p[k+6] .* xx
        effects[4, :, i] = p[k+7] .* xx
        effects[5, :, i] .= p[k+8]
        k+=9
    end
    effects[6, :, :] .= floor(p[37]) #nG1
    effects[7, :, :] .= floor(p[38]) #nG2
    effects[8, :, :] .= floor(p[39]) #nD1
    effects[9, :, :] .= floor(p[40]) #nD2

    return effects
end

function residHillAll(hillParams::Array{Float64,1}, concentrations::Array{Float64,2}, g1::Array{Float64,3}, g2::Array{Float64,3})
    res = Atomic{eltype(hillParams)}(0.0)
    params = getODEparamsAll(hillParams, concentrations)

    # Solve for all drugs
    @threads for j = 1:4
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
function optimize_hillAll(concs::Array{Float64,2}, g1::Array{Float64,3}, g2::Array{Float64,3}; maxstep = 1E5)
    hillCostAll(hillParams) = residHillAll(hillParams, concs, g1, g2)

    low = [minimum(concs[:, 1]), 0.01, 1e-9, 1e-9, 1e-9, 1e-9, 0.0, 0.0, 0.45, minimum(concs[:, 2]), 0.01, 1e-9, 1e-9, 1e-9, 1e-9, 0.0, 0.0, 0.45, minimum(concs[:, 3]), 0.01, 1e-9, 1e-9, 1e-9, 1e-9, 0.0, 0.0, 0.45, minimum(concs[:, 4]), 0.01, 1e-9, 1e-9, 1e-9, 1e-9, 0.0, 0.0, 0.45, 2, 10, 0, 0]
    high = [maximum(concs[:, 1]), 10.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 0.55, maximum(concs[:, 2]), 10.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 0.55, maximum(concs[:, 3]), 10.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 0.55, maximum(concs[:, 4]), 10.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 0.55, 60, 180, 50, 50]

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