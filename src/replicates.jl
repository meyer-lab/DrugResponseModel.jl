""" This file includes functions to calculate values related to replicates. """

### A function to predict G1 and G2 for the three replicates.###
function predict_replicates(p1, p2, p3, g0)
    t = LinRange(0.0, 95.0, 189)
    ps = cat(p1, p2, p3, dims = 3)
    G1 = ones(189, 8, 3)
    G2 = ones(189, 8, 3)

    for i = 1:size(ps, 2) # concentration number
        for j = 1:size(ps, 3)
            G1[:, i, j], G2[:, i, j], _ = predict(ps[:, i, j], g0, t)
        end
    end

    return G1, G2 # all simulation
end


""" A function to calculate std and mean of ODE parameters for each drug. """
function mean_std_params(effs1, effs2, effs3)
    eff = cat(effs1, effs2, effs3, dims = 3)
    return mean(eff, dims = 3), std(eff, dims = 3)
end


""" plots all the three simulated trials. Along with avg and std of data. """
function plot2(G1, G2, g1s1, g1s2, g1s3, g2s1, g2s2, g2s3, conc::Float64, i::Int, j::Int)
    time = LinRange(0.0, 95.0, 189)
    G1s = cat(g1s1, g1s2, g1s3, dims = 4)
    G2s = cat(g2s1, g2s2, g2s3, dims = 4)
    meang1 = mean(G1s, dims = 4)
    meang2 = mean(G2s, dims = 4)
    stdg1 = std(G1s, dims = 4)
    stdg2 = std(G2s, dims = 4)

    Plots.plot(
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
    Plots.plot!(time, G1[:, 1], label = "G1", color = 6)
    Plots.plot!(time, G1[:, 2], label = "", color = 6)
    Plots.plot!(time, G1[:, 3], label = "", color = 6)
    Plots.plot!(time, meang2[:, i, j]; ribbon = stdg2[:, i, j], color = 7, label = "", alpha = 0.1)
    Plots.plot!(time, G2[:, 1], label = "G2", color = 7)
    Plots.plot!(time, G2[:, 2], label = "", color = 7)
    Plots.plot!(time, G2[:, 3], label = "", color = 7)
    ylims!((0.0, 45))
end
