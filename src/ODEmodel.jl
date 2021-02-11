"""
    In this file we want to estimate parameters of an ODE model describing the
    number of cells in G1 or G2 phase of the cell cycle 
"""

const nG1 = 8
const nG2 = 20
const nSp = nG1 + nG2

""" Make the transition matrix. """
function ODEjac(p::AbstractVector{T}, t::Real)::Matrix{T} where {T <: Real}
    # with 2 G1 and 4 G2 we have p = [a1, a2, b1, b2, b3, b4, g11, g12, g21, g22, g23, g24] # 12 params
    #                         indx = [1,  2,  3,  4,  5,  6,  7,   8,   9,   10,  11,  12 ]
    A = zeros(nSp, nSp)

    A[diagind(A, 0)[1:Int(nG1/2)]] .= -(p[1] + p[7])
    A[diagind(A, 0)[Int(nG1/2 + 1):nG1]] .= -(p[2] + p[8])
    A[diagind(A, 0)[(nG1 + 1):Int(nG1 + nG2/4)]] .= -(p[3] + p[9])
    A[diagind(A, 0)[Int(nG1 + nG2/4 + 1):Int(nG1 + nG2/2)]] .= -(p[4] + p[10])
    A[diagind(A, 0)[Int(nG1 + nG2/2 + 1):Int(nG1 + 3*nG2/4)]] .= -(p[5] + p[11])
    A[diagind(A, 0)[Int(nG1 + 3*nG2/4 + 1):nSp]] .= -(p[6] + p[12])

    A[diagind(A, -1)[1:Int(nG1/2)]] .= p[1]
    A[diagind(A, -1)[Int(nG1/2 + 1):nG1]] .= p[2]
    A[diagind(A, -1)[(nG1 + 1):Int(nG1 + nG2/4)]] .= p[3]
    A[diagind(A, -1)[Int(nG1 + nG2/4 + 1):Int(nG1 + nG2/2)]] .= p[4]
    A[diagind(A, -1)[Int(nG1 + nG2/2 + 1):Int(nG1 + 3*nG2/4)]] .= p[5]
    A[diagind(A, -1)[Int(nG1 + 3*nG2/4 + 1):(nSp-1)]] .= p[6]

    A[1, nSp] = 2 * p[6]

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

    return 20.0 * v / sum(v)
end


function vTOg(v::AbstractVector)
    G1 = sum(view(v, 1:nG1))
    G2 = sum(view(v, (nG1 + 1):nSp))
    return G1, G2
end


""" Predicts the model given a set of parametrs. """
function predict(p::AbstractVector, g_0::AbstractVector, t::Union{Real, LinRange}, g1data = nothing, g2data = nothing)
    @assert length(p) == 12 # we have 2 G1 prog rates, 4 G2 prog rates, 2 G1 death and 4 G2 death rates.

    if length(g_0) == length(p)
        v = startV(g_0)
    else
        @assert length(g_0) == nSp
        v = copy(g_0)
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
