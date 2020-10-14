""" This file includes functions to calculate values related to replicates. """

### A function to predict G1 and G2 for the three replicates.###
function predict_replicates(p1, p2, p3, concs1, g0)

    t = LinRange(0.0, 95.0, 189)
    G1_1 = ones(189, 8)
    G2_1 = ones(189, 8)
    G1_2 = ones(189, 8)
    G2_2 = ones(189, 8)
    G1_3 = ones(189, 8)
    G2_3 = ones(189, 8)

    for i = 1:8 # concentration number
        G1_1[:, i], G2_1[:, i], _ = predict(p1[:, i], g0, t)
        G1_2[:, i], G2_2[:, i], _ = predict(p2[:, i], g0, t)
        G1_3[:, i], G2_3[:, i], _ = predict(p3[:, i], g0, t)
    end

    return G1_1, G2_1, G1_2, G2_2, G1_3, G2_3 # all simulation
end

""" A function to calculate std and mean of ODE parameters for each drug. """
function mean_std_params(effs1, effs2, effs3)
    meann = ones(9, 8)
    stdd = ones(9, 8)
    for i = 1:8
        for j = 1:9
            meann[j, i] = mean([effs1[j, i], effs2[j, i], effs3[j, i]])
            stdd[j, i] = std([effs1[j, i], effs2[j, i], effs3[j, i]])
        end
    end
    return meann, stdd
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


""" A Function to find mean and std of data in G1 and G2 separately. """
function mean_std_data(G1_1, G1_2, G1_3, G2_1, G2_2, G2_3)

    meanG1 = ones(189, 8)
    meanG2 = ones(189, 8)
    stdG1 = ones(189, 8)
    stdG2 = ones(189, 8)
    for j = 1:8
        for k = 1:189
            meanG1[k, j] = mean([G1_1[k, j], G1_2[k, j], G1_3[k, j]])
            meanG2[k, j] = mean([G2_1[k, j], G2_2[k, j], G2_3[k, j]])
            stdG1[k, j] = std([G1_1[k, j], G1_2[k, j], G1_3[k, j]])
            stdG2[k, j] = std([G2_1[k, j], G2_2[k, j], G2_3[k, j]])
        end
    end
    return meanG1, meanG2, stdG1, stdG2
end

""" plots all the three simulated trials. Along with avg and std of data. """
function plot2(G1_1, G1_2, G1_3, G2_1, G2_2, G2_3, g1s1, g1s2, g1s3, g2s1, g2s2, g2s3, conc, i, j)
    time = LinRange(0.0, 95.0, 189)
    meang1 = zeros(189, 8, 5)
    meang2 = zeros(189, 8, 5)
    stdg1 = zeros(189, 8, 5)
    stdg2 = zeros(189, 8, 5)
    for k = 1:5
        meang1[:, :, k], meang2[:, :, k], stdg1[:, :, k], stdg2[:, :, k] =
            mean_std_data(g1s1[:, :, k], g1s2[:, :, k], g1s3[:, :, k], g2s1[:, :, k], g2s2[:, :, k], g2s3[:, :, k])
    end
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
