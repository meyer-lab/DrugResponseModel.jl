"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle 
"""

""" Make the transition matrix. """
function ODEjac(p::Vector{T}, nG1::Int, nG2::Int, nD1::Int, nD2::Int) where {T}
    # p = [alpha, beta, gamma1, gamma2, nG1, nG2, nD1, nD2]
    v1 = [-ones(nG1) * (p[3] + p[1]); -ones(nG2) * (p[4] + p[2])]
    v2 = [ones(nG1) * p[1]; ones(nG2 - 1) * p[2]]

    if nD1 == 1
        append!(v1, -ones(nD1) * p[3])
        append!(v2, 0.0)
    elseif nD1 > 1
        append!(v1, -ones(nD1) * p[3])
        append!(v2, [0.0; ones(nD1 - 1) * p[3]])
    end

    if nD2 == 1
        append!(v1, -ones(nD2) * p[4])
        append!(v2, 0.0)
    elseif nD2 > 1
        append!(v1, -ones(nD2) * p[4])
        append!(v2, [0.0; ones(nD2 - 1) * p[4]])
    end

    A = diagm(0 => v1, -1 => v2)

    A[1, nG1 + nG2] = 2 * p[2]
    if nD1 > 0
        A[nG1 + nG2 + 1, 1:nG1] = p[3] * ones(1, nG1)
    end
    if nD2 > 0
        A[nG1 + nG2 + nD1 + 1, (nG1 + 1):(nG1 + nG2)] = p[4] * ones(1, nG2)
    end

    return A
end


function vTOg(v::Vector, nG1::Int, nG2::Int, nD1::Int, nD2::Int)
    G1 = sum(view(v, 1:nG1)) + sum(view(v, (nG1 + nG2 + 1):(nG1 + nG2 + nD1)))
    G2 = sum(view(v, (nG1 + 1):(nG1 + nG2))) + sum(view(v, (nG1 + nG2 + nD1 + 1):(nG1 + nG2 + nD1 + nD2)))
    return G1, G2
end


""" Predicts the model given a set of parametrs. """
function predict(p, g_0, t, g1data = nothing, g2data = nothing)
    # Convert parameters to phase numbers
    nG1, nG2, nD1, nD2 = Tuple(Int.(floor.(p[6:9])))

    if g_0 isa Real
        v = [ones(nG1) * p[5] * g_0 / nG1; ones(nG2) * (1.0 - p[5]) * g_0 / nG2]

        if nD1 + nD2 > 0
            append!(v, zeros(nD1 + nD2))
        end
    else
        v = g_0
    end

    A = ODEjac(p, nG1, nG2, nD1, nD2)

    if t isa Real
        v = ExponentialUtilities.expv(t, A, v)

        G1, G2 = vTOg(v, nG1, nG2, nD1, nD2)
    else
        # Some assumptions
        @assert t[1] == 0.0
        rmul!(A, t[2])

        if eltype(p) === Float64
            A = LinearAlgebra.exp!(A)
        else
            A = ExponentialUtilities.exp_generic(A)
        end

        if g1data === nothing
            G1 = Vector{eltype(p)}(undef, length(t))
            G2 = Vector{eltype(p)}(undef, length(t))

            for ii = 1:length(t)
                G1[ii], G2[ii] = vTOg(v, nG1, nG2, nD1, nD2)

                v = A * v
            end
        else
            cost = 0.0

            for ii = 1:length(t)
                G1, G2 = vTOg(v, nG1, nG2, nD1, nD2)
                cost += norm(G1 - g1data[ii]) + norm(G2 - g2data[ii])

                v = A * v
            end

            return cost, v
        end
    end

    return G1, G2, v
end
