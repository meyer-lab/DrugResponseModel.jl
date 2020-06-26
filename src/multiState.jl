""" Make the transition matrix. """

function ODEjacMulti(p::Vector{T}, nG1::Int, nG2::Int, nD1::Int, nD2::Int) where {T}
    # p = [alpha1, alpha2, alpha3, alpha4, beta1, beta2, gamma1, gamma2, gamma3, gamma4, nG1, nG2, nD1, nD2]
    # nGi is the same across the states, so it nDi.
    # the order is:
    # nG1 --> nG2 --> nG1 --> nG2
    if nD1 == 0
        D1 = T[]
        diagD1 = T[]
    elseif nD1 == 1
        D1 = [0.0]
        diagD1 = -ones(nD1) * p[3]
    else
        D1 = [0.0; ones(nD1 - 1) * p[3]]
        diagD1 = -ones(nD1) * p[3]
    end

    if nD2 == 0
        D2 = T[]
        diagD2 = T[]
    elseif nD2 == 1
        D2 = [0.0]
        diagD2 = -ones(nD2) * p[4]
    else
        D2 = [0.0; ones(nD2 - 1) * p[4]]
        diagD2 = -ones(nD2) * p[4]
    end

    v1 = [-ones(nG1) * (p[7] + p[1]); -ones(nG2) * (p[8] + p[2]); -ones(nG1) * (p[9] + p[3]); -ones(nG2) * (p[10] + p[4]); diagD1; diagD2; diagD1; diagD2]
    v2 = [ones(nG1) * p[1]; ones(nG2 - 1) * p[2]; ones(nG1) * p[3]; ones(nG2 - 1) * p[4]; D1; D2; D1; D2]
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
