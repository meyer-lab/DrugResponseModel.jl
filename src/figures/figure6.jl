""" plot hill curves for drug effects. """

function plot_hillCurves()

    concs, _, _, _ = load(189, 1)
    p = parameters()
    ps = getODEparams(p, concs)

    # average parameters
    new_PS = zeros(4, 8, 5)
    new_PS[1, :, :] .= mean(ps[1:2, :, :], dims=1)[1, :, :]
    new_PS[2, :, :] .= mean(ps[3:6, :, :], dims=1)[1, :, :]
    new_PS[3, :, :] .= mean(ps[7:8, :, :], dims=1)[1, :, :]
    new_PS[4, :, :] .= mean(ps[9:12, :, :], dims=1)[1, :, :]

    p0 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p1=plot(concs[:, 1], new_PS[:, :, 1]', labels=["G1 prog." "G2 prog." "G1 death" "G2 death"], title="lapatinib effects", xlabel="concentration[nM]",lw=2, ylabel="effect", titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.4cm,
    fg_legend = :transparent,
    top_margin = 1.4cm,
    left_margin = 1.3cm,
    right_margin = 1.3cm)
    p2=plot(concs[:, 2], new_PS[:, :, 2]', labels=["G1 prog." "G2 prog." "G1 death" "G2 death"], title="doxorubicin effects", xlabel="concentration[nM]",lw=2, ylabel="effect", titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.4cm,
    fg_legend = :transparent,
    top_margin = 1.4cm,
    left_margin = 1.3cm,
    right_margin = 1.3cm,)
    p3=plot(concs[:, 3], new_PS[:, :, 3]', labels=["G1 prog." "G2 prog." "G1 death" "G2 death"], title="gemcitabine effects", xlabel="concentration[nM]",lw=2, ylabel="effect", titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.4cm,
    fg_legend = :transparent,
    top_margin = 1.4cm,
    left_margin = 1.3cm,
    right_margin = 1.3cm,)
    p4=plot(concs[:, 4], new_PS[:, :, 4]', labels=["G1 prog." "G2 prog." "G1 death" "G2 death"], title="paclitaxel effects", xlabel="concentration[nM]",lw=2, ylabel="effect", titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.4cm,
    fg_legend = :transparent,
    top_margin = 1.4cm,
    left_margin = 1.3cm,
    right_margin = 1.3cm,)
    p5=plot(concs[:, 5], new_PS[:, :, 5]', labels=["G1 prog." "G2 prog." "G1 death" "G2 death"], title="palbociclib effects", xlabel="concentration[nM]", lw=2, ylabel="effect", titlefont = Plots.font("Helvetica", 12),
    legendfont = Plots.font("Helvetica", 9),
    guidefont = Plots.font("Helvetica", 12),
    xtickfont = Plots.font("Helvetica", 12),
    ytickfont = Plots.font("Helvetica", 12),
    bottom_margin = 1.4cm,
    fg_legend = :transparent,
    top_margin = 1.4cm,
    left_margin = 1.3cm,
    right_margin = 1.3cm,)

    p=plot(p0, p1, p2, p3, p4, p5, size=(1200, 700), layout=(2, 3))
    savefig(p, "figure6.svg")
end