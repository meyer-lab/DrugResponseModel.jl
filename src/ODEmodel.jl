"""
    In this file we want to estimate parameters of an ODE model describing the
    number of cells in G1 or G2 phase of the cell cycle 
"""

const nG1 = 6
const nG2 = 6
const nSp = nG1 + nG2

""" Make the transition matrix. """
function ODEjac(p::AbstractVector{T}, t::Real) where {T <: Real}
    # p = [alpha1, alpha2, beta1, beta2, gamma11, gamma12, gamma21, gamma22, %G1]
    v1 = [-ones(Int(nG1/2)) * (p[5] + p[1]); -ones(Int(nG1/2)) * (p[2] + p[6]); -ones(Int(nG2/2)) * (p[3] + p[7]); -ones(Int(nG2/2)) * (p[4] + p[8])]
    v2 = [ones(Int(nG1/2)) * p[1]; ones(Int(nG1/2)) * p[2]; ones(Int(nG2/2)) * p[3]; ones(Int(nG2/2 - 1)) * p[4]]

    A = diagm(0 => v1, -1 => v2)
    A[1, nSp] = 2 * p[4]

    lmul!(t, A)
    return LinearAlgebra.exp!(A)
end


""" Find the starting vector from the steady-state of the control condition. """
function startV(p::AbstractVector{T})::AbstractVector{T} where {T <: Real}
    _, _, v = predict(p, 1.0, 100.0)

    for ii = 1:100
        v /= sum(v)
        _, _, v = predict(p, v, 100.0)
    end

    return v / sum(v)
end


function vTOg(v::AbstractVector)
    G1 = sum(view(v, 1:nG1))
    G2 = sum(view(v, (nG1 + 1):nSp))
    return G1, G2
end

""" To find the initial vector. """
function residue(p, g1, g2)
    t = LinRange(0.0, 0.5 * size(g1, 1), size(g1, 1))
    g00 = g1[1, 1] + g2[1, 1]
    return predict(p, g00, t, g1[:, 1], g2[:, 1])[1]
end

function optimise_control(g1, g2)
    smin = [1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 0.25];
    smax = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.75];
    f(x) = residue(x, g1, g2)
    return DrugResponseModel.optimize_helper(f, smin, smax, 100000)
end

function newPredict(p, pControl, t::Union{Real, LinRange}, g1data = nothing, g2data = nothing)

    vStart = startV(pControl)
    # Note that vStart is always scaled so the starting cell number is 1.0

    return predict(p, vStart, t, g1data, g2data)
end


""" Predicts the model given a set of parametrs. """
function predict(p::AbstractVector, g_0, t::Union{Real, LinRange}, g1data = nothing, g2data = nothing)
    @assert length(p) == 9
    @assert all(p .>= 0.0)

    if g_0 isa Real
        v = SVector{nSp}([ones(nG1) * p[9] * g_0 / nG1; ones(nG2) * (1.0 - p[9]) * g_0 / nG2])
    else
        v = g_0
    end

    if t isa Real
        A = ODEjac(p, t)

        v = A * v
        G1, G2 = vTOg(v)
    else
        # Some assumptions
        @assert t.start == 0.0
        A = ODEjac(p, t[2])

        if g1data === nothing
            G1 = Vector{eltype(p)}(undef, length(t))
            G2 = Vector{eltype(p)}(undef, length(t))

            for ii = 1:length(t)
                G1[ii], G2[ii] = vTOg(v)

                v = A * v
            end
        else
            cost = 0.0

            for ii = 1:length(t)
                G1, G2 = vTOg(v)
                cost += norm(G1 - g1data[ii]) + norm(G2 - g2data[ii])

                v = A * v
            end

            return cost, v
        end
    end

    return G1, G2, v
end
