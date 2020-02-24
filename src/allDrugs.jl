""" In this file we fit all the drugs att once. """

function getODEparamsAll(p::Vector, concentrations::Matrix)
    effects = Matrix{eltype(p)}(undef, 9, 8, 4)

    k = 1
    # Scaled drug effect
    for j=1:4
        xx = 1.0 ./ (1.0 .+ (p[k] ./ concentrations[:, j]) .^ p[k+1])

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

function residHillAll(hillParams::Vector, concentrations::Matrix, g1::Matrix, g2::Matrix)
    res = Atomic{eltype(hillParams)}(0.0)
    params = getODEparamsAll(hillParams, concentrations)

    # Solve for all drugs
    @threads for j = 1:4
        @threads for ii = 1:length(concentrations[:, j])
            atomic_add!(
                res,
                cost(
                    params[:, ii, j],
                    g1[:, ii],
                    g2[:, ii],
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