""" Make the transition matrix. """

function ODEjacMulti(p::Vector{T}, nG1::Int, nG2::Int, nD1::Int, nD2::Int) where {T}
    # p = [alpha1, alpha2, alpha3, alpha4, beta1, beta2, gamma1, gamma2, gamma3, gamma4, nG1, nG2, nD1, nD2]
    # nGi is the same across the states, so it nDi.
    # the order is:
    # nG1 --> nG2 --> nG1 --> nG2 after that:
    # nD1 --> nD2 --> nD1 --> nD2

    if nD1 == 0
        D1 = T[]
        diagD1 = T[]
        D3 = T[]
        diagD3 = T[]
    elseif nD1 == 1
        D1 = [0.0]
        diagD1 = -ones(nD1) * p[7]
        D3 = [0.0]
        diagD3 = -ones(nD1) * p[9]
    else
        D1 = [0.0; ones(nD1 - 1) * p[7]]
        diagD1 = -ones(nD1) * p[7]
        D3 = [0.0; ones(nD1 - 1) * p[9]]
        diagD3 = -ones(nD1) * p[9]
    end

    # nD2
    if nD2 == 0
        D2 = T[]
        diagD2 = T[]
        D4 = T[]
        diagD4 = T[]
    elseif nD2 == 1
        D2 = [0.0]
        diagD2 = -ones(nD2) * p[8]
        D4 = [0.0]
        diagD4 = -ones(nD2) * p[10]
    else
        D2 = [0.0; ones(nD2 - 1) * p[8]]
        diagD2 = -ones(nD2) * p[8]
        D4 = [0.0; ones(nD2 - 1) * p[10]]
        diagD4 = -ones(nD2) * p[10]
    end

    v1 = [-ones(nG1) * (p[7] + p[1]); -ones(nG2) * (p[8] + p[2]); -ones(nG1) * (p[9] + p[3]); -ones(nG2) * (p[10] + p[4]); diagD1; diagD2; diagD3; diagD4]
    v2 = [ones(nG1) * p[1]; ones(nG2 - 1) * p[2]; ones(nG1) * p[3]; ones(nG2 - 1) * p[4]; D1; D2; D3; D4]
    A = diagm(0 => v1, -1 => v2)

    A[1, nG1 + nG2] = 2 * p[2]
    A[1, 2*nG1 + 2*nG2] = 2 * p[6]
    A[nG1 + nG2, nG1 + nG2] -= p[5]
    A[2*nG1 + 2*nG2, 2*nG1 + 2*nG2] -= p[6]
    A[nG1 + nG2 + 1, 2*nG1 + 2*nG2] = 2 * p[4]
    A[2*nG1 + 2*nG2, 2*nG1 + nG2 + 1] = 2 * p[4]
    if nD1 > 0
        A[2*nG1 + 2*nG2 + 1, 1:nG1] = p[7] * ones(1, nG1)
        A[2*nG1 + 2*nG2 + nD1 + nD2 + 1, (nG1 + nG2 + 1):(2*nG1 + nG2)] = p[9] * ones(1, nG1)
    end
    if nD2 > 0
        A[2*nG1 + 2*nG2 + nD1 + 1, (nG1 + 1):(nG1 + nG2)] = p[8] * ones(1, nG2)
        A[2*nG1 + 2*nG2 + 2*nD1 + nD2 + 1, (2*nG1 + nG2 + 1):(2*nG1 + 2*nG2)] = p[10] * ones(1, nG2)
    end

    return A
end

""" Predicts the model given a set of parametrs. """
function Multipredict(p, g_0, t)
    # Convert parameters to phase numbers
    nG1 = Int(floor(p[11]))
    nG2 = Int(floor(p[12]))
    nD1 = Int(floor(p[13]))
    nD2 = Int(floor(p[14]))
    if nD1 == 0
        D1 = Float64[]
        D3 = Float64[]
    else
        D1 = zeros(nD1)
        D3 = zeros(nD1)
    end
    if nD2 == 0
        D2 = Float64[]
        D4 = Float64[]
    else
        D2 = zeros(nD2)
        D4 = zeros(nD2)
    end

    if g_0 isa Real
        v = [ones(nG1) * p[15] * g_0 / nG1; ones(nG2) * (1.0 - p[15]) * g_0 / nG2; ones(nG1) * p[15] * g_0 / nG1; ones(nG2) * (1.0 - p[15]) * g_0 / nG2; D1; D2; D3; D4]
    else
        v = g_0
    end

    A = ODEjacMulti(p, nG1, nG2, nD1, nD2)

    if t isa Real
        v = ExponentialUtilities.expv(t, A, v)

        G1 = sum(v[1:nG1]) + sum(v[(nG1 + nG2 + 1):(nG1 + nG2 + nD1)])
        G2 = sum(v[(nG1 + 1):(nG1 + nG2)]) + sum(v[(nG1 + nG2 + nD1 + 1):(nG1 + nG2 + nD1 + nD2)])
    else
        # Some assumptions
        @assert t[1] == 0.0
        rmul!(A, t[2])
        A = LinearAlgebra.exp!(A)

        G11 = Vector{eltype(p)}(undef, length(t))
        G12 = Vector{eltype(p)}(undef, length(t))
        G21 = Vector{eltype(p)}(undef, length(t))
        G22 = Vector{eltype(p)}(undef, length(t))

        for ii = 1:length(G11)
            G11[ii] = sum(view(v, 1:nG1)) + sum(view(v, (2*nG1 + 2*nG2 + 1):(2*nG1 + 2*nG2 + nD1)))
            G12[ii] = sum(view(v, (nG1 + 1):(nG1 + nG2))) + sum(view(v, (2*nG1 + 2*nG2 + nD1 + 1):(2*nG1 + 2*nG2 + nD1 + nD2)))
            G21[ii] = sum(view(v, (nG1 + nG2):(2*nG1 + nG2))) + sum(view(v, (2*nG1 + 2*nG2 + nD1 + nD2 + 1):(2*nG1 + 2*nG2 + 2*nD1 + nD2)))
            G22[ii] = sum(view(v, (2*nG1 + nG2 + 1):(2*nG1 + 2*nG2))) + sum(view(v, (2*nG1 + 2*nG2 + 2*nD1 + nD2 + 1):(2*nG1 + 2*nG2 + 2*nD1 + 2*nD2)))

            v = A * v
        end
    end

    return G11, G12, G21, G22, v
end


""" Calculates the cost function for a given set of parameters. """
function Multicost(p, g11, g12, g21, g22)
    t = LinRange(0.0, 0.5 * length(g1), length(g1))
    G11, G12, G21, G22 = Multipredict(p, g11[1] + g12[1] + g21[1] + g22[1], t)

    return norm(G11 - g11) + norm(G12 - g12) + norm(G21 - g21) + norm(G22 - g22)
end
