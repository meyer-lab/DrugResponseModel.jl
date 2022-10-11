""" Plot cell line fits when fitting ALL AT ONCE."""

""" Helper for plotting time series fits. """
function plot_fig1(concs, g1, g1data, tite, G, subPlabel, palet, time)
    p = Plots.plot(
        time,
        g1,
        lw = 4,
        legend = :topleft,
        label = ["control" "$(concs[4]) nM" "$(concs[5]) nM" "$(concs[6]) nM" "$(concs[7]) nM" "$(concs[8]) nM"],
        fg_legend = :transparent,
        palette = palet,
        title = tite,
        titlefont = Plots.font("Helvetica", 14),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 14),
        xtickfont = Plots.font("Helvetica", 14),
        ytickfont = Plots.font("Helvetica", 14),
        xlabel = "time [hr]",
        xticks = 0:24.0:96.0,
        ylabel = "$G cell number",
        bottom_margin = 1.25cm,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    Plots.plot!(time, g1data, lw = 4, linestyle = :dot, label = ["" "" "" "" "" "" ""])
    annotate!(-0.5, 1.5, Plots.text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 4.0))
    p
end

""" Helper function to plot accumulated dead cells. """
function plot_dead_acc(concs, drugs, efcs, siz)

    g = zeros(siz, 8, 4, 8) # total
    t = LinRange(0.0, 96.0, siz)
    d = zeros(siz, 6, 4) # datapoints x concs x drugs
    ls = [1, 4, 5, 6, 7, 8]
    for i = 1:4
        k = 1
        for j in ls
            g[:, j, i, 1], g[:, j, i, 2], g[:, j, i, 3], g[:, j, i, 4], g[:, j, i, 5], g[:, j, i, 6], g[:, j, i, 7], g[:, j, i, 8] =
                DrugResponseModel.predictD(efcs[:, j, i], efcs[:, 1, i], t)
            d[:, k, i] =
                efcs[9, j, i] .* g[:, j, i, 1] .+ efcs[10, j, i] .* g[:, j, i, 2] .+ efcs[11, j, i] .* g[:, j, i, 3] .+
                efcs[12, j, i] .* g[:, j, i, 4] .+ efcs[13, j, i] .* g[:, j, i, 5] .+ efcs[14, j, i] .* g[:, j, i, 6] .+
                efcs[15, j, i] .* g[:, j, i, 7] .+ efcs[16, j, i] .* g[:, j, i, 8]
            k += 1
        end
    end
    intg = zeros(siz, 6, 4)
    for i = 1:4
        for j = 1:6
            intg[:, j, i] = cumul_integrate(t, d[:, j, i])
        end
    end
    function single_dead(intt, i)
        p1 = Plots.plot(
            t,
            intt[:, :, i],
            labels =  ["control" "$(concs[i][4]) nM" "$(concs[i][5]) nM" "$(concs[i][6]) nM" "$(concs[i][7]) nM" "$(concs[i][8]) nM"],
            title = drugs[i],
            lw = 2,
            ylabel = "Accumulated cell death #",
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
        ylims!((-0.05, 3.0))
        return p1
    end

    p = [single_dead(intg, j) for j=1:4]
    return p

end

ps_hcc = [1.36845, 4.05584, 3.99669, 3.99714, 0.102503, 0.458295, 0.659276, 3.99775, 0.365036, 0.627148, 6.94891e-5, 0.000136151, 0.0127684, 0.0517801, 3.20979e-5, 4.19056e-5, 7.37436e-6, 0.0228549, 1.11308, 0.264787, 3.99759, 1.45882e-5, 3.06579e-6, 0.375592, 0.102024, 3.99633, 0.320449, 0.722804, 3.53551e-6, 0.032066, 1.67137e-7, 4.27832e-7, 8.83771e-8, 2.89182e-6, 4.5214e-7, 6.8609e-7, 16.4934, 2.45814, 1.12299, 1.55116, 0.0033046, 0.185649, 0.673943, 0.305188, 0.0594824, 0.0924852, 5.03661e-5, 7.25378e-6, 0.00764243, 3.2136e-6, 1.2219e-6, 7.56255e-7, 0.00440004, 0.00037535, 0.741965, 1.89296, 3.98903, 3.99547, 2.31225e-6, 0.383961, 0.464131, 2.02298, 0.091986, 0.293161, 1.52185e-5, 2.38363e-5, 0.021285, 1.42916e-6, 3.24183e-8, 3.92282e-6, 1.48194e-7, 1.03137e-6, 3.99914, 3.99967, 0.482352, 0.358652, 0.531748, 3.99842, 0.35508, 0.537322]
ps_mt1 = [1.57556, 1.62607, 4.35979e-6, 2.20023, 6.40848e-9, 1.63139e-7, 4.0, 1.03588, 2.82295e-7, 1.518, 7.26753e-9, 3.32831e-8, 7.38826e-9, 0.0847927, 1.59243e-8, 6.84661e-9, 0.0635307, 2.0618e-9, 18.4637, 0.498753, 3.99999, 4.0, 0.00330672, 0.509792, 4.0, 0.669307, 3.92732, 1.30048, 9.84433e-8, 1.17613e-8, 1.03609e-9, 1.63366e-8, 4.43533e-8, 1.11375e-9, 1.15564e-8, 2.01859e-9, 41.7844, 1.70309, 2.01363e-6, 0.0480213, 1.40979e-9, 1.5706e-7, 3.99999, 0.89413, 4.0, 7.42191e-9, 1.00757e-6, 0.0640029, 7.28897e-9, 1.06527, 2.77222e-8, 4.83221e-9, 4.92092e-9, 1.18316e-9, 0.733685, 1.82624, 4.0, 0.0108699, 0.991755, 0.0497405, 3.99663, 0.614765, 2.77109e-9, 1.10172, 1.66923e-7, 0.00238703, 6.67847e-8, 1.01807e-9, 5.22232e-9, 2.70425e-9, 1.80354e-9, 2.04378e-9, 4.0, 4.0, 0.574277, 2.19519, 3.9953, 0.627565, 4.0, 1.23897]
ps_mda = [2.13839, 2.83884, 1.51441, 1.50356, 0.0863139, 0.289618, 0.753405, 7.01123e-6, 0.0869785, 3.99967, 1.53553e-5, 7.3733e-6, 1.41833e-7, 0.0528917, 0.0221832, 0.374224, 2.10944e-7, 7.27666e-8, 412.122, 0.685044, 5.93963e-6, 1.44106e-6, 0.269538, 0.0777158, 0.345028, 3.99898, 0.305618, 3.99972, 0.0880853, 0.0556922, 5.90658e-7, 6.21931e-7, 6.23376e-8, 3.69143e-7, 3.19237e-8, 4.25434e-6, 12.8393, 1.10132, 4.33747e-7, 3.75999e-7, 0.258907, 4.06523e-6, 3.42362e-7, 0.190433, 0.267539, 0.18205, 0.0253506, 0.0219874, 0.00864397, 0.00222323, 0.00597668, 4.48199e-7, 1.28088e-6, 1.45763e-7, 0.471157, 5.8535, 2.15688, 1.73648, 0.253272, 0.13052, 0.160907, 0.129681, 0.464879, 3.99994, 0.0109992, 0.0343582, 9.64949e-7, 9.13167e-8, 1.38366e-8, 0.01639, 0.0109108, 0.000553012, 3.99994, 3.99988, 0.104727, 1.34732, 0.475374, 3.99997, 0.298439, 3.99996]

""" To plot the fits and accumulated cell death for each cell line, we do the following:
1. tensor, names, concs, conds = DrugResponseModel.__cellLineName__()
where __cellLineName__ could be one of [hcc_all, mt1_all, mda_all]
2. imporing the estimated parameters according to the cell line, one of [ps_hcc, ps_mt1, ps_mda] above.
3. DrugResponseModel.figure70(tensor, names, concs, conds, ps)"""

function figure70(tensor, names, concs, conds, ps)
    ENV["GKSwstype"]="nul"
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end
    efcs = getODEparams(ps, cs)
    t = LinRange(0.0, 96, size(tensor)[2])
    Gshort = zeros(size(tensor)[2], 6, 4, 2)
    cc = [1, 4, 5, 6, 7, 8]

    for (i, drug) in enumerate(names)
        for j=1:6
            Gshort[:, j, i, 1], Gshort[:, j, i, 2], _ = DrugResponseModel.predict(efcs[:, cc[j], i], efcs[:, 1, i], t)
        end
    end

    # params at EC50
    ec50 = zeros(16, 4)
    conc_ec50 = zeros((1, 4))
    k=1
    for i=1:4
        conc_ec50[1, i] = ps[k]
        k+=18
    end
    ec50 = getODEparams(ps, conc_ec50)[:, 1, :]

    # phase durations
    # @ control
    gi = zeros(2, 2, 4) # g1/g2 x control/ec50 x drugs
    gi[1, 1, :] .= (2 ./ efcs[1, 1, :] .+ 2 ./ efcs[2, 1, :] .+ 2 ./ efcs[3, 1, :] .+ 2 ./ efcs[4, 1, :])
    gi[2, 1, :] .= (5 ./ efcs[5, 1, :] .+ 5 ./ efcs[6, 1, :] .+ 5 ./ efcs[7, 1, :] .+ 5 ./ efcs[8, 1, :])

    # @ ec50
    gi[1, 2, :] .= (2 ./ ec50[1, :] .+ 2 ./ ec50[2, :] .+ 2 ./ ec50[3, :] .+ 2 ./ ec50[4, :])
    gi[2, 2, :] .= (5 ./ ec50[5, :] .+ 5 ./ ec50[6, :] .+ 5 ./ ec50[7, :] .+ 5 ./ ec50[8, :])

    gmshort = zeros(size(tensor)[2], 6, 4, 2) # datapoints x concs x drugs x g1/g2
    for i=1:2
        gmshort[:, 1, :, i] .= tensor[i, :, 1, :]
        gmshort[:, 2:6, :, i] .= tensor[i, :, 4:8, :]
    end

    # cell deaths
    deathContG1 = zeros(4, 4)
    deathEC50G1 = zeros(4, 4)
    deathContG1[1, :] .= (efcs[9, 1, :]) ./ (efcs[9, 1, :] .+ efcs[1, 1, :])
    deathEC50G1[1, :] .= (ec50[9, :]) ./ (ec50[9, :] .+ ec50[1, :])
    for i = 2:4
        deathContG1[i, :] .= (1 .- deathContG1[i - 1, :]) .* (efcs[i + 8, 1, :]) ./ (efcs[i + 8, 1, :] .+ efcs[i, 1, :])
        deathEC50G1[i, :] .= (1 .- deathEC50G1[i - 1, :]) .* (ec50[i + 8, :]) ./ (ec50[i + 8, :] .+ ec50[i, :])
    end
    deathContG2 = zeros(4, 4)
    deathEC50G2 = zeros(4, 4)
    deathContG2[1, :] .= (efcs[13, 1, :]) ./ (efcs[13, 1, :] .+ efcs[5, 1, :])
    deathEC50G2[1, :] .= (ec50[13, :]) ./ (ec50[13, :] .+ ec50[5, :])
    for i = 14:16
        deathContG2[i - 12, :] = (1 .- deathContG2[i - 13, :]) .* (efcs[i, 1, :]) ./ (efcs[i, 1, :] .+ efcs[i - 8, 1, :])
        deathEC50G2[i - 12, :] = (1 .- deathEC50G2[i - 13, :]) .* (ec50[i, :]) ./ (ec50[i, :] .+ ec50[i - 8, :])
    end

    p0 = Plots.plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p1 = plot_fig1(concs[1], Gshort[:, :, 1, 1], gmshort[:, :, 1, 1], string("Dynamical Model Fits - ", names[1]), "G1", "A", :PuBu_6, t)
    p2 = plot_fig1(concs[1], Gshort[:, :, 1, 2], gmshort[:, :, 1, 2], string("Dynamical Model Fits - ", names[1]), "S-G2", "B", :PuBu_6, t)
    p3 = plot_fig1(concs[2], Gshort[:, :, 2, 1], gmshort[:, :, 2, 1], string("Dynamical Model Fits - ", names[2]), "G1", "C", :PuBu_6, t)
    p4 = plot_fig1(concs[2], Gshort[:, :, 2, 2], gmshort[:, :, 2, 2], string("Dynamical Model Fits - ", names[2]), "S-G2", "D", :PuBu_6, t)
    p5 = plot_fig1(concs[3], Gshort[:, :, 3, 1], gmshort[:, :, 3, 1], string("Dynamical Model Fits - ", names[3]), "G1", "E", :PuBu_6, t)
    p6 = plot_fig1(concs[3], Gshort[:, :, 3, 2], gmshort[:, :, 3, 2], string("Dynamical Model Fits - ", names[3]), "S-G2", "F", :PuBu_6, t)
    p7 = plot_fig1(concs[4], Gshort[:, :, 4, 1], gmshort[:, :, 4, 1], string("Dynamical Model Fits - ", names[4]), "G1", "G", :PuBu_6, t)
    p8 = plot_fig1(concs[4], Gshort[:, :, 4, 2], gmshort[:, :, 4, 2], string("Dynamical Model Fits - ", names[4]), "S-G2", "H", :PuBu_6, t)
    p = plot_dead_acc(concs, names, efcs, size(Gshort)[1])

    figure1 = Plots.plot(p1, p2, p[1], p3, p4, p[2], p5, p6, p[3], p7, p8, p[4], size = (1500, 1800), layout = (4, 3))
    Plots.savefig(figure1, "figure70.svg")
end
