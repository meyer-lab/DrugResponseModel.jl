"""Calculate deadcells in each subphase."""

function vTOgD(v::AbstractVector)
    G11 = sum(view(v, 1:2))
    G12 = sum(view(v, 3:4))
    G13 = sum(view(v, 5:6))
    G14 = sum(view(v, 7:8))
    G21 = sum(view(v, 9:13))
    G22 = sum(view(v, 14:18))
    G23 = sum(view(v, 19:23))
    G24 = sum(view(v, 24:28))
    return G11, G12, G13, G14, G21, G22, G23, G24
end

""" Predicts the model given a set of parametrs. """
function predictD(p::AbstractVector, g_0::AbstractVector, t::Union{Real, LinRange}, g1data = nothing, g2data = nothing)
    @assert length(p) == 16 # we have 2 G1 prog rates, 4 G2 prog rates, 2 G1 death and 4 G2 death rates.

    if length(g_0) == length(p)
        v = startV(g_0)
    else
        @assert length(g_0) == nSp
        v = copy(g_0)
    end

    # Some assumptions
    @assert t.start == 0.0
    A = ODEjac(p)
    lmul!(t[2], A)
    A = LinearAlgebra.exp!(A)
    u = similar(v)

    G11 = Vector{eltype(p)}(undef, length(t))
    G12 = Vector{eltype(p)}(undef, length(t))
    G13 = Vector{eltype(p)}(undef, length(t))
    G14 = Vector{eltype(p)}(undef, length(t))
    G21 = Vector{eltype(p)}(undef, length(t))
    G22 = Vector{eltype(p)}(undef, length(t))
    G23 = Vector{eltype(p)}(undef, length(t))
    G24 = Vector{eltype(p)}(undef, length(t))

    for ii = 1:length(t)
        G11[ii], G12[ii], G13[ii], G14[ii], G21[ii], G22[ii], G23[ii], G24[ii] = vTOgD(v)

        mul!(u, A, v)
        copyto!(v, u)
    end

    return G11, G12, G13, G14, G21, G22, G23, G24
end

function output_deadcells()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    ps = parameters()
    efcs = getODEparams(ps, concs)
    g = zeros(189, 8, 5, 8) # total
    t = LinRange(0.0, 96.0, 189)
    d = zeros(189, 6, 5)
    ls = [1, 4, 5, 6, 7, 8]
    for i = 1:5
        k = 1
        for j in ls
            g[:, j, i, 1], g[:, j, i, 2], g[:, j, i, 3], g[:, j, i, 4], g[:, j, i, 5], g[:, j, i, 6], g[:, j, i, 7], g[:, j, i, 8] =
                predictD(efcs[:, j, i], efcs[:, 1, i], t)
            d[:, k, i] =
                efcs[9, j, i] .* g[:, j, i, 1] .+ efcs[10, j, i] .* g[:, j, i, 2] .+ efcs[11, j, i] .* g[:, j, i, 3] .+
                efcs[12, j, i] .* g[:, j, i, 4] .+ efcs[13, j, i] .* g[:, j, i, 5] .+ efcs[14, j, i] .* g[:, j, i, 6] .+
                efcs[15, j, i] .* g[:, j, i, 7] .+ efcs[16, j, i] .* g[:, j, i, 8]
            k += 1
        end
    end
    intg = zeros(189, 6, 5)
    for i = 1:5
        for j = 1:6
            intg[:, j, i] = cumul_integrate(t, d[:, j, i])
        end
    end
    p1 = Plots.plot(
        t,
        intg[:, :, 1],
        labels = ["control" "25nM" "50nM" "100nM" "250nM" "500nM"],
        title = "lapatinib",
        lw = 2,
        ylabel = "accumulated cell death #",
        xlabel = "time [hr]",
        palette = :YlOrRd_6,
        legend = :left,
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    ylims!((-0.05, 2.0))

    p3 = Plots.plot(
        t,
        intg[:, :, 3],
        labels = ["control" "2.5nM" "5nM" "10nM" "30nM" "100nM"],
        title = "gemcitabine",
        lw = 2,
        ylabel = "accumulated cell death #",
        xlabel = "time [hr]",
        palette = :YlOrRd_6,
        legend = :left,
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.5cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    ylims!((-0.05, 2.0))

    return p1, p3
end
