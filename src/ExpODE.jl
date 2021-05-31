"""
    In this file we want to estimate parameters of an ODE model describing the
    number of cells in G1 or G2 phase of the cell cycle 
"""

const nG1 = 1
const nG2 = 1
const nSp = nG1 + nG2

""" Make the transition matrix. """
function ODEjac_exp(p::AbstractVector{T})::Matrix{T} where {T <: Real}
    # with 2 G1 and 4 G2 we have p = [a1, a2, b1, b2, b3, b4, g11, g12, g21, g22, g23, g24] # 12 params
    #                        index = [1,  2,  3,  4,  5,  6,  7,   8,   9,   10,  11,  12 ]
    #                                [a, 0, b, 0, 0, 0, g1, 0, g2, 0, 0, 0]
    A = zeros(2, 2)
    A[1] = -(p[1] + p[7])
    A[2] = 2 * p[3]
    A[3] = p[1]
    A[4] = -(p[3] + p[9])

    return A
end

""" Find the starting vector from the steady-state of the control condition. """
function startV_exp(p::AbstractVector{T})::AbstractVector{T} where {T <: Real}
    @assert all(p .>= 0.0)
    @assert all(p[7:end] .== 0.0) # No cell death in the control
    A = ODEjac_exp(p)

    vals, vecs = eigen(A)

    selectt = real.(vals) .> 0.0
    @assert sum(selectt) == 1
    vecs = vec(vecs[:, selectt])
    @assert all(isreal.(vals[selectt]))
    @assert all(isreal.(vecs))

    return vecs / sum(vecs)
end


function vTOg_exp(v::AbstractVector)
    G1 = sum(view(v, 1:nG1))
    G2 = sum(view(v, (nG1 + 1):nSp))
    return G1, G2
end


""" Predicts the model given a set of parametrs. """
function predict_exp(p::AbstractVector, g_0::AbstractVector, t::Union{Real, LinRange}, g1data = nothing, g2data = nothing)
    @assert length(p) == 12 # we have 2 G1 prog rates, 4 G2 prog rates, 2 G1 death and 4 G2 death rates.

    if length(g_0) == length(p)
        v = startV_exp(g_0)
    else
        @assert length(g_0) == nSp
        v = copy(g_0)
    end

    if t isa Real
        A = ODEjac_exp(p)
        lmul!(t, A)
        A = LinearAlgebra.exp!(A)

        v = A * v
        G1, G2 = vTOg_exp(v)
    else
        # Some assumptions
        @assert t.start == 0.0
        A = ODEjac_exp(p)
        lmul!(t[2], A)
        A = LinearAlgebra.exp!(A)
        u = similar(v)

        if g1data === nothing
            G1 = Vector{eltype(p)}(undef, length(t))
            G2 = Vector{eltype(p)}(undef, length(t))

            for ii = 1:length(t)
                G1[ii], G2[ii] = vTOg_exp(v)

                mul!(u, A, v)
                copyto!(v, u)
            end
        else
            cost = 0.0

            for ii = 1:length(t)
                G1, G2 = vTOg_exp(v)
                cost += norm(G1 - g1data[ii]) + norm(G2 - g2data[ii])

                mul!(u, A, v)
                copyto!(v, u)
            end

            return cost, v
        end
    end

    return G1, G2, v
end

function residHillexp(x::Vector, conc::Vector, g1::Matrix, g2::Matrix)

    res = 0.0
    for i=3:8
        res += 60*(maximum([0, (x[i] - x[i + 12])]))^2
    end
    params = getODEparams(x, conc)
    t = LinRange(0.0, 0.5 * size(g1, 1), size(g1, 1))
    # Solve each concentration separately
    for ii = 1:length(conc)
        res += predict_exp(params[:, ii, 1], params[:, 1, 1], t, g1[:, ii], g2[:, ii])[1]
    end
    return res
end

function residHillAllexp(hP, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0

    # Solve for all drugs
    t = 1
    for j = 1:5
        hill = hP[[t:(t + 13); 71:76]]
        res += residHillexp(hill, concentrations[:, j], g1[:, :, j], g2[:, :, j])
        t += 14
    end

    return res
end

function optim_allexp(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 800000)
    f(x) = residHillAllexp(x, concs, g1, g2)

    # [a, 0, b, 0, 0, 0, g1, 0, g2, 0, 0, 0]
    lP = [minimum(concs); 0.01; 1e-9; 0; 1e-9; 0.0; 0.0; 0.0; 1e-9; 0; 1e-9; 0.0; 0.0; 0.0]
    low = vcat(lP, lP, lP, lP, lP, 1e-9, 0, 1e-9, 0, 0, 0)
    hP = [maximum(concs); 10.0; 2; 0.0; 2; 0.0; 0.0; 0; 2; 0.0; 2; 0.0; 0.0; 0.0]
    high = vcat(hP, hP, hP, hP, hP, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0)

    return optimize_helper(f, low, high, maxiter)
end
# p_optim = [69.0898, 1.75704, 1.00525e-9, 0.0, 1.30328e-9, 0.0, 0.0, 0.0, 1.01896e-9, 0.0, 0.0117861, 0.0, 0.0, 0.0, 52.0114, 1.04672, 0.0249361, 0.0, 1.20034e-9, 0.0, 0.0, 0.0, 1.0814e-9, 0.0, 0.0336839, 0.0, 0.0, 0.0, 6.48332, 1.5805, 0.232475, 0.0, 0.0762636, 0.0, 0.0, 0.0, 0.205384, 0.0, 0.00899069, 0.0, 0.0, 0.0, 3.3574, 3.51098, 0.0213074, 0.0, 1.03062e-9, 0.0, 0.0, 0.0, 1.00861e-9, 0.0, 0.0212547, 0.0, 0.0, 0.0, 29.4276, 1.27817, 0.334571, 0.0, 4.21098e-9, 0.0, 0.0, 0.0, 1.01008e-9, 0.0, 0.00287931, 0.0, 0.0, 0.0, 0.206861, 0.0, 0.0145702, 0.0, 0.0, 0.0]