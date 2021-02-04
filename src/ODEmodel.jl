"""
    In this file we want to estimate parameters of an ODE model describing the
    number of cells in G1 or G2 phase of the cell cycle 
"""

const nG1 = 5
const nG2 = 5
const nSp = nG1 + nG2

""" Make the transition matrix. """
function ODEjac(p::Vector{T}, t::Real) where {T}
    # p = [alpha, beta, gamma1, gamma2]
    v1 = [-ones(nG1) * (p[3] + p[1]); -ones(nG2) * (p[4] + p[2])]
    v2 = [ones(nG1) * p[1]; ones(nG2 - 1) * p[2]]

    A = diagm(0 => v1, -1 => v2)
    A[1, nSp] = 2 * p[2]

    A = SMatrix{nSp, nSp}(A)
    A = exp(A * t)

    return A
end


function vTOg(v::AbstractVector)
    G1 = sum(view(v, 1:nG1))
    G2 = sum(view(v, (nG1 + 1):nSp))
    return G1, G2
end


""" Predicts the model given a set of parametrs. """
function predict(p, g_0, t::Union{Real, LinRange}, g1data = nothing, g2data = nothing)

    if g_0 isa Real
        v = SVector{nSp}([ones(nG1) * p[5] * g_0 / nG1; ones(nG2) * (1.0 - p[5]) * g_0 / nG2])
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
