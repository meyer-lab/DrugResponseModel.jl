"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle 
"""
global nG1 = 8;
global nG2 = 24;

""" Make the transition matrix. """
function ODEjac(p::Vector{T}) where {T}
    # p = [alpha1, alpha2, beta1, beta2, gamma11, gamma12, gamma21, gamma22]
    v1 = [-ones(Int(nG1/2)) * (p[5] + p[1]); -ones(Int(nG1/2)) * (p[2] + p[6]); -ones(Int(nG2/2)) * (p[3] + p[7]); -ones(Int(nG2/2)) * (p[4] + p[8])]
    v2 = [ones(Int(nG1/2)) * p[1]; ones(Int(nG1/2)) * p[2]; ones(Int(nG2/2 - 1)) * p[3]; ones(Int(nG2/2 - 1)) * p[4]]

    A = diagm(0 => v1, -1 => v2)

    A[1, nG1 + nG2] = 2 * p[4]
    A[nG1 + nG2, end-1] = p[4]

    return A
end

function vTOg(v::Vector, nG1::Int, nG2::Int)
    G1 = sum(view(v, 1:nG1))
    G2 = sum(view(v, (nG1 + 1):(nG1 + nG2)))
    return G1, G2
end


""" Predicts the model given a set of parametrs. """
function predict(p, g_0, t, g1data = nothing, g2data = nothing)
    # Convert parameters to phase numbers

    if g_0 isa Real
        v = [ones(nG1) * p[9] * g_0 / nG1; ones(nG2) * (1.0 - p[9]) * g_0 / nG2]

    else
        v = g_0
    end

    A = ODEjac(p[1:8])

    if t isa Real
        rmul!(A, t)
        A = LinearAlgebra.exp!(A)

        v = A * v
        G1, G2 = vTOg(v, nG1, nG2)
    else
        # Some assumptions
        @assert t[1] == 0.0
        rmul!(A, t[2])

        A = LinearAlgebra.exp!(A)

        if g1data === nothing
            G1 = Vector{eltype(p)}(undef, length(t))
            G2 = Vector{eltype(p)}(undef, length(t))

            for ii = 1:length(t)
                G1[ii], G2[ii] = vTOg(v, nG1, nG2)

                v = A * v
            end
        else
            cost = 0.0

            for ii = 1:length(t)
                G1, G2 = vTOg(v, nG1, nG2)
                cost += norm(G1 - g1data[ii]) + norm(G2 - g2data[ii])

                v = A * v
            end

            return cost, v
        end
    end

    return G1, G2, v
end
