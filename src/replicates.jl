""" This file includes functions to calculate values related to replicates. """

### A function to predict G1 and G2 for the three replicates.###
function predict_replicates(p1, p2, p3, g0)
    t = LinRange(0.0, 95.0, 189)
    ps = cat(p1, p2, p3, dims = 3)
    G1 = ones(189, 8, 3)
    G2 = ones(189, 8, 3)

    for i = 1:size(ps, 2) # concentration number
        for j = 1:size(ps, 3)
            G1[:, i, j], G2[:, i, j], _ = predict(p1[:, i, j], g0, t)
        end
    end

    return G1[:, :, 1], G1[:, :, 2], G1[:, :, 3], G2[:, :, 1], G2[:, :, 2], G2[:, :, 3] # all simulation
end

""" A function to calculate std and mean of ODE parameters for each drug. """
function mean_std_params(effs1, effs2, effs3)
    eff = cat(effs1, effs2, effs3, dims = 3)
    return mean(eff, dims = 3), std(eff, dims = 3)
end

""" Calculate the # of cells in G1 for a set of parameters and T """
function numcells(params, g0)
    @assert(all(params .>= 0.0), "negative params $params")
    t = LinRange(0.0, 94.5, 189)
    G1, G2 = predict(params, g0, t)

    @assert(all(G1[2:end] .>= 0.0), "negative cell number in G1 $G1")
    @assert(all(G2[2:end] .>= 0.0), "negative cell number in G2 $G2")
    return G1[end] + G2[end]
end


""" plots all the three simulated trials. Along with avg and std of data. """
function plot2(G1_1, G1_2, G1_3, G2_1, G2_2, G2_3, g1s1, g1s2, g1s3, g2s1, g2s2, g2s3, conc, i, j)
    time = LinRange(0.0, 95.0, 189)
    G1 = cat(g1s1, g1s2, g1s3, dims = 4)
    G2 = cat(g2s1, g2s2, g2s3, dims = 4)
    meang1 = mean(G1, dims = 4)
    meang2 = mean(G2, dims = 4)
    stdg1 = std(G1, dims = 4)
    stdg2 = std(G2, dims = 4)

    plot(
        time,
        meang1[:, i, j];
        ribbon = stdg1[:, i, j],
        title = string(conc, "nM"),
        color = 6,
        label = "",
        xlabel = "time [hr]",
        ylabel = "cell number",
        alpha = 0.1,
    )
    plot!(time, G1_1, label = "G1", color = 6)
    plot!(time, G1_2, label = "", color = 6)
    plot!(time, G1_3, label = "", color = 6)
    plot!(time, meang2[:, i, j]; ribbon = stdg2[:, i, j], color = 7, label = "", alpha = 0.1)
    plot!(time, G2_1, label = "G2", color = 7)
    plot!(time, G2_2, label = "", color = 7)
    plot!(time, G2_3, label = "", color = 7)
    ylims!((0.0, 45))
end
