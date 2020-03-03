""" In this file we fit all the drugs att once. """

function getODEparamsAll(p::Array{Float64, 1}, concentrations::Array{Float64, 2})
    effects = effects = zeros(9, 8, 4)

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
    effects[5, :, :] .= p[27]
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

    low = [
        minimum(concs[:, 1]),
        0.01,
        1e-9,
        1e-9,
        0.0,
        0.0,
        minimum(concs[:, 2]),
        0.01,
        1e-9,
        1e-9,
        0.0,
        0.0,
        minimum(concs[:, 3]),
        0.01,
        1e-9,
        1e-9,
        0.0,
        0.0,
        minimum(concs[:, 4]),
        0.01,
        1e-9,
        1e-9,
        0.0,
        0.0,
        1e-9,
        1e-9,
        0.45,
        2,
        10,
        0,
        0,
    ]
    high = [
        maximum(concs[:, 1]),
        10.0,
        3.0,
        3.0,
        1.0,
        1.0,
        maximum(concs[:, 2]),
        10.0,
        3.0,
        3.0,
        1.0,
        1.0,
        maximum(concs[:, 3]),
        10.0,
        3.0,
        3.0,
        1.0,
        1.0,
        maximum(concs[:, 4]),
        10.0,
        3.0,
        3.0,
        1.0,
        1.0,
        3.0,
        3.0,
        0.55,
        60,
        180,
        50,
        50,
    ]

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


""" Combination functions. """
function ParamForBliss(p)
    """ To calculate Bliss independence drug effect
    we assume delays are constant, death rates are additive,
    and will keep the alpha and beta intact."""
    par = zeros(3,8)
    par[1,:] = p[1,:] # alpha stays the same
    par[2,:] = p[2,:] # beta stays the same
    par[3,:] = p[3,:] + p[4,:] # additivity assumption for death rates
    @assert(all(par .>= 0.0), "You cannot use Bliss combination, the individual drug effects <= 0")
    @assert(all(par .<= 1.0), "You cannot use Bliss combination, the individual drug effects >= 1.0")
    return par
end

function BlissCombination(p1::Matrix{Float64}, p2::Matrix{Float64})
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect independently. """

    param1 = ParamForBliss(p1)
    param2 = ParamForBliss(p2)
    """ For 8x8 combination of drug concentrations, and the 3 is for G1 progression rate, G2 progression rate, and death rates. """
    deathBliss = zeros(8,8,3)
    for j in 1:8
        for k in 1:8
            deathBliss[j,k,:] .= param1[:,j] .+ param2[:,k] .- param1[:,j] .* param2[:,k]
            end
        end
    return param1, param2, deathBliss
end