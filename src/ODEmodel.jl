"""
    In this file we want to estimate parameters of an ODE model describing the
    number of cells in G1 or G2 phase of the cell cycle 
"""

const nG1 = 6
const nG2 = 6
const nSp = nG1 + nG2

""" Make the transition matrix. """
function ODEjac(p::AbstractVector{T}, t::Real)::Matrix{T} where {T <: Real}
    # p = [alpha, beta, gamma1, gamma2]
    A = zeros(nSp, nSp)

    A[diagind(A, 0)[1:nG1]] .= -(p[3] + p[1])
    A[diagind(A, 0)[(nG1 + 1):end]] .= -(p[4] + p[2])

    A[diagind(A, -1)[1:nG1]] .= p[1]
    A[diagind(A, -1)[(nG1 + 1):end]] .= p[2]

    A[1, nSp] = 2 * p[2]

    lmul!(t, A)
    return LinearAlgebra.exp!(A)
end


""" Find the starting vector from the steady-state of the control condition. """
function startV(p::AbstractVector{T})::AbstractVector{T} where {T <: Real}
    v = ones(nSp)
    A = ODEjac(p, 100.0)
    u = similar(v)

    for ii = 1:100
        v /= sum(v)
        mul!(u, A, v)
        copyto!(v, u)
    end

    return v / sum(v)
end


function vTOg(v::AbstractVector)
    G1 = sum(view(v, 1:nG1))
    G2 = sum(view(v, (nG1 + 1):nSp))
    return G1, G2
end


""" Predicts the model given a set of parametrs. """
function predict(p, g_0, t::Union{Real, LinRange}, g1data = nothing, g2data = nothing)
    if g_0 isa Real
        v = Vector{eltype(p)}(undef, nSp)
        v[1:nG1] .= p[5] * g_0 / nG1
        v[(nG1 + 1):end] .= (1.0 - p[5]) * g_0 / nG2
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
        u = similar(v)

        if g1data === nothing
            G1 = Vector{eltype(p)}(undef, length(t))
            G2 = Vector{eltype(p)}(undef, length(t))

            for ii = 1:length(t)
                G1[ii], G2[ii] = vTOg(v)

                mul!(u, A, v)
                copyto!(v, u)
            end
        else
            cost = 0.0

            for ii = 1:length(t)
                G1, G2 = vTOg(v)
                cost += norm(G1 - g1data[ii]) + norm(G2 - g2data[ii])

                mul!(u, A, v)
                copyto!(v, u)
            end

            return cost, v
        end
    end

    return G1, G2, v
end
