""" plot hill curves for drug effects. """

function plot_hillCurves()

    concs, _, _, _ = load(189, 1)
    p = parameters()
    ps = getODEparams(p, concs)

    # average parameters
    durations = zeros(2, 8, 5)
    durations[1, :, :] .= (2 ./ ps[1, :, :] .+ 2 ./ ps[2, :, :] .+ 2 ./ ps[3, :, :] .+ 2 ./ ps[4, :, :])
    durations[2, :, :] .= (5 ./ ps[5, :, :] .+ 5 ./ ps[6, :, :] .+ 5 ./ ps[7, :, :] .+ 5 ./ ps[8, :, :])
    deathG1 = zeros(4, 8, 5)
    deathG1[1, :, :] .= (ps[9, :, :]) ./ (ps[9, :, :] .+ ps[1, :, :])
    for i = 2:4
        deathG1[i, :, :] .= (1 .- deathG1[i - 1, :, :]) .* (ps[i + 8, :, :]) ./ (ps[i + 8, :, :] .+ ps[i, :, :])
    end
    deathG2 = zeros(4, 8, 5)
    deathG2[1, :, :] .= (ps[13, :, :]) ./ (ps[13, :, :] .+ ps[5, :, :])
    for i = 14:16
        deathG2[i - 12, :, :] = (1 .- deathG2[i - 13, :, :]) .* (ps[i, :, :]) ./ (ps[i, :, :] .+ ps[i - 8, :, :])
    end
    death = zeros(2, 8, 5)

    for i = 1:5
        death[1, :, i] = sum(deathG1[:, :, i], dims=1)
        death[2, :, i] = sum(deathG2[:, :, i], dims=1)
    end

    function plot_efs(concs, new_PS, drugname, ylabel, ec50, y)

        plot(
        concs,
        new_PS',
        labels = ["G1" "G2"],
        title = "$drugname effects",
        xlabel = "Concentration[nM]",
        lw = 3,
        ylabel = ylabel,
        titlefont = Plots.font("Helvetica", 15),
        legendfont = Plots.font("Helvetica", 11),
        guidefont = Plots.font("Helvetica", 15),
        xtickfont = Plots.font("Helvetica", 15),
        ytickfont = Plots.font("Helvetica", 15),
        bottom_margin = 1.4cm,
        fg_legend = :transparent,
        top_margin = 1.4cm,
        left_margin = 1.3cm,
        right_margin = 1.3cm,
        xrotation=30,
    )
    plot!([ec50], seriestype = :vline, alpha=0.6, label="EC50", lw=3)
    ylims!((0.0, y))
    end

    p0 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p1 = plot_efs(concs[:, 1], durations[:, :, 1], "Lapatinib", "Avg. phase duration", p[1], 300.0)
    p2 = plot_efs(concs[:, 2], durations[:, :, 2], "Doxorubicin", "Avg. phase duration", p[19], 300.0)
    p3 = plot_efs(concs[:, 3], durations[:, :, 3], "Gemcitabine", "Avg. phase duration", p[37], 300.0)
    p4 = plot_efs(concs[:, 4], durations[:, :, 4], "Paclitaxel", "Avg. phase duration", p[55], 300.0)
    p5 = plot_efs(concs[:, 5], durations[:, :, 5], "Palbociclib", "Avg. phase duration", p[73], 300.0)
    p6 = plot_efs(concs[:, 1], death[:, :, 1], "Lapatinib", "Avg. cell death prob.", p[1], 0.6)
    p7 = plot_efs(concs[:, 2], death[:, :, 2], "Doxorubicin", "Avg. cell death prob.", p[19], 0.6)
    p8 = plot_efs(concs[:, 3], death[:, :, 3], "Gemcitabine", "Avg. cell death prob.", p[37], 0.6)
    p9 = plot_efs(concs[:, 4], death[:, :, 4], "Paclitaxel", "Avg. cell death prob.", p[55], 0.6)
    p10 = plot_efs(concs[:, 5], death[:, :, 5], "Palbociclib", "Avg. cell death prob.", p[73], 0.6)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, size = (1800, 700), layout = (2, 5))
    savefig(p, "figure6.svg")
end
