""" Plot hcc1143 fits when fitting ALL AT ONCE."""
function figure7()
    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    params = [6.76402, 1.45485, 0.0901376, 2.00322, 2.00283, 0.187222, 2.005, 2.00412, 0.148397, 0.0586825, 0.0213053, 2.52195e-7, 2.41545e-8, 1.44707e-8, 3.29222e-8, 1.37208e-8, 5.59402e-9, 0.00126715, 43.3336, 1.39287, 0.139558, 0.138035, 0.0904822, 0.0493381, 9.62971e-8, 5.98781e-7, 0.264581, 0.00759821, 6.58923e-9, 2.90942e-8, 7.08595e-8, 4.87951e-8, 2.21321e-8, 1.84993e-8, 0.0515671, 6.56555e-9, 1.05139, 2.02354, 0.0172485, 2.00418, 0.258644, 0.949032, 2.002, 2.00531, 1.1072e-7, 0.119229, 2.07555e-8, 0.407274, 8.36152e-8, 9.54029e-8, 2.32011e-7, 4.88497e-8, 0.0021182, 4.3538e-8, 1.14876, 3.86715, 0.15917, 2.01371, 2.01491, 1.44785, 2.00591, 2.01811, 0.140134, 0.749386, 0.0433281, 1.6227e-7, 5.84897e-8, 5.88594e-8, 2.78218e-8, 1.63851e-8, 5.49249e-9, 6.74174e-8, 5000.0, 0.310183, 0.172035, 2.00422, 2.00337, 1.4227, 2.00219, 2.00051, 0.0141024, 0.236535, 0.0479571, 2.37249e-7, 1.46444e-7, 5.86096e-8, 2.4064e-8, 1.12806e-8, 7.34672e-9, 8.34734e-8, 35.3323, 1.03572, 0.162309, 2.00718, 2.00352, 1.41681, 2.00072, 2.0004, 0.0604097, 0.662226, 1.67197e-8, 1.02143e-7, 5.22464e-8, 2.72402e-8, 2.40392e-8, 7.06357e-9, 8.59539e-9, 0.0646765, 0.0750058, 2.0, 1.99998, 1.41318, 2.0, 2.0, 0.126361, 0.810005]
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end
    pode = getODEparams(params, cs)

    pp = []
    for (i, drug) in enumerate(names)

        Gs = zeros(182, 8, 2)
        t = LinRange(0.0, 95, 182)

        for j=1:8
            Gs[:, j, 1], Gs[:, j, 2], _ = DrugResponseModel.predict(pode[:, j, i], pode[:, 1, i], t)
        end

        push!(pp, DrugResponseModel.Eachdrug_sim(Gs[:, :, 1], tensor[1, :, :, i], concs[i], "G1", drug))
        push!(pp, DrugResponseModel.Eachdrug_sim(Gs[:, :, 2], tensor[2, :, :, i], concs[i], "S/G2", drug))
    end

    pl = plotGrid((3, 4), [pp...];)
    return draw(SVG("figure7.svg", 18inch, 14inch), pl)

end

function plot_fig1(concs, g1, g1data, tite, G, subPlabel, palet, time)

    p = Plots.plot(
        time,
        g1,
        lw = 2,
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
    Plots.plot!(time, g1data, lw = 2, linestyle = :dot, label = ["" "" "" "" "" "" ""])
    annotate!(-1.0, 2.0, Plots.text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 2.5))
    p
end

function plot_pG1_mean(y1, y2, ymax, Phasename, ylabel, subPlabel, plus)

    x = ["Lapatinib", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib"]
    conts = Plots.scatter(
        x,
        y1',
        color = "red",
        xlabel = "drugs",
        xrotation = 30,
        label = "Control",
        markerstrokewidth = 0,
        markersize = 10,
        ylabel = ylabel,
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
        title = "$Phasename effects",
    )
    Plots.scatter!(
        x,
        y2',
        color = "cyan4",
        xlabel = "drugs",
        label = "EC50",
        markerstrokewidth = 0,
        markersize = 10,
        ylabel = ylabel,
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
        title = "$Phasename effects",
    )
    annotate!(-0.5, (ymax + plus), Plots.text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((-0.1, ymax))
end


function figure70()
    ENV["GKSwstype"]="nul"
    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    ps = [6.76402, 1.45485, 0.0901376, 2.00322, 2.00283, 0.187222, 2.005, 2.00412, 0.148397, 0.0586825, 0.0213053, 2.52195e-7, 2.41545e-8, 1.44707e-8, 3.29222e-8, 1.37208e-8, 5.59402e-9, 0.00126715, 43.3336, 1.39287, 0.139558, 0.138035, 0.0904822, 0.0493381, 9.62971e-8, 5.98781e-7, 0.264581, 0.00759821, 6.58923e-9, 2.90942e-8, 7.08595e-8, 4.87951e-8, 2.21321e-8, 1.84993e-8, 0.0515671, 6.56555e-9, 1.05139, 2.02354, 0.0172485, 2.00418, 0.258644, 0.949032, 2.002, 2.00531, 1.1072e-7, 0.119229, 2.07555e-8, 0.407274, 8.36152e-8, 9.54029e-8, 2.32011e-7, 4.88497e-8, 0.0021182, 4.3538e-8, 1.14876, 3.86715, 0.15917, 2.01371, 2.01491, 1.44785, 2.00591, 2.01811, 0.140134, 0.749386, 0.0433281, 1.6227e-7, 5.84897e-8, 5.88594e-8, 2.78218e-8, 1.63851e-8, 5.49249e-9, 6.74174e-8, 5000.0, 0.310183, 0.172035, 2.00422, 2.00337, 1.4227, 2.00219, 2.00051, 0.0141024, 0.236535, 0.0479571, 2.37249e-7, 1.46444e-7, 5.86096e-8, 2.4064e-8, 1.12806e-8, 7.34672e-9, 8.34734e-8, 35.3323, 1.03572, 0.162309, 2.00718, 2.00352, 1.41681, 2.00072, 2.0004, 0.0604097, 0.662226, 1.67197e-8, 1.02143e-7, 5.22464e-8, 2.72402e-8, 2.40392e-8, 7.06357e-9, 8.59539e-9, 0.0646765, 0.0750058, 2.0, 1.99998, 1.41318, 2.0, 2.0, 0.126361, 0.810005]
    cs = zeros(8, length(concs))
    for i=1:length(concs)
        cs[:, i] = concs[i]
    end
    efcs = getODEparams(ps, cs)
    t = LinRange(0.0, 95, 182)
    Gshort = zeros(182, 6, 6, 2)
    cc = [1, 4, 5, 6, 7, 8]

    for (i, drug) in enumerate(names)
        for j=1:6
            Gshort[:, j, i, 1], Gshort[:, j, i, 2], _ = DrugResponseModel.predict(efcs[:, cc[j], i], efcs[:, 1, i], t)
        end
    end

    # params at EC50
    ec50 = zeros(16, 6)
    conc_ec50 = zeros((1, 6))
    k=1
    for i=1:6
        conc_ec50[1, i] = ps[k]
        k+=18
    end
    ec50 = getODEparams(ps, conc_ec50)[:, 1, :]

    # phase durations
    # @ control
    gi = zeros(2, 2, 6) # g1/g2 x control/ec50 x drugs
    gi[1, 1, :] .= (2 ./ efcs[1, 1, :] .+ 2 ./ efcs[2, 1, :] .+ 2 ./ efcs[3, 1, :] .+ 2 ./ efcs[4, 1, :])
    gi[2, 1, :] .= (5 ./ efcs[5, 1, :] .+ 5 ./ efcs[6, 1, :] .+ 5 ./ efcs[7, 1, :] .+ 5 ./ efcs[8, 1, :])

    # @ ec50
    gi[1, 2, :] .= (2 ./ ec50[1, :] .+ 2 ./ ec50[2, :] .+ 2 ./ ec50[3, :] .+ 2 ./ ec50[4, :])
    gi[2, 2, :] .= (5 ./ ec50[5, :] .+ 5 ./ ec50[6, :] .+ 5 ./ ec50[7, :] .+ 5 ./ ec50[8, :])

    gmshort = zeros(182, 6, 6, 2) # datapoints x concs x drugs x g1/g2
    for i=1:2
        gmshort[:, 1, :, i] .= tensor[i, :, 1, :]
        gmshort[:, 2:6, :, i] .= tensor[i, :, 4:8, :]
    end

    p0 = Plots.plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p1 = plot_fig1(concs[1], Gshort[:, :, 1, 1], gmshort[:, :, 1, 1], "Dynamical Model Fits - BEZ235", "G1", "A", :PuBu_6, t)
    p2 = plot_fig1(concs[1], Gshort[:, :, 1, 2], gmshort[:, :, 1, 2], "Dynamical Model Fits - BEZ235", "S/G2", "B", :PuBu_6, t)
    p3 = plot_fig1(concs[2], Gshort[:, :, 2, 1], gmshort[:, :, 2, 1], "Dynamical Model Fits - Doxorubicin", "G1", "C", :PuBu_6, t)
    p4 = plot_fig1(concs[2], Gshort[:, :, 2, 2], gmshort[:, :, 2, 2], "Dynamical Model Fits - Doxorubicin", "S/G2", "D", :PuBu_6, t)
    p5 = plot_pG1_mean(gi[1, 1, :]', gi[1, 2, :]', 60.0, "G1 ", "Avg phase duration [hr]", "E", 4.0)
    p6 = plot_pG1_mean(gi[2, 1, :]', gi[2, 2, :]', 60.0, "S/G2 ", "Avg phase duration [hr]", "F", 4.0)
    # p7 = plot_pG1_mean(sum(deathContG1, dims = 1), sum(deathEC50G1, dims = 1), 0.2, "G1", "Cell death probability", "G", 1.0)
    # p8 = plot_pG1_mean(sum(deathContG2, dims = 1), sum(deathEC50G2, dims = 1), 0.2, "S/G2", "Cell death probability", "H", 1.0)


    figure1 = Plots.plot(p1, p2, p3, p4, p5, p6, p0,p0, size = (2000, 1000), layout = (2, 4))
    Plots.savefig(figure1, "figure70.svg")
end
