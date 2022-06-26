""" Plot hcc1143 fits when fitting ALL AT ONCE."""
function figure7()
    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    params = [9.20945, 2.12023, 2.98291, 0.430111, 0.09537, 0.0956342, 0.242984, 0.0828915, 1.71385, 0.299598, 0.000505937, 0.0647114, 5.01559e-6, 1.49885e-5, 0.00544226, 4.35115e-6, 3.03577e-5, 2.15171e-5, 13.9863, 2.77415, 3.06754, 3.09629, 0.129791, 0.002272, 2.12184, 0.461891, 0.362844, 0.0293403, 0.000990191, 0.00063532, 2.98906e-5, 3.40495e-6, 0.000103298, 5.23445e-7, 7.58079e-6, 0.00361932, 1.11528, 1.93197, 3.82999, 3.8369, 7.05901e-6, 3.52973e-5, 0.714193, 0.0650093, 0.521932, 0.0465728, 0.00111718, 6.65118e-6, 0.00254141, 0.0378996, 9.51571e-5, 7.83845e-7, 0.000242097, 0.00102866, 1.5658, 3.63417, 3.88274, 3.97669, 0.0956175, 0.64807, 2.13597, 0.396049, 1.0315, 0.395203, 0.00162117, 0.00053917, 0.000125216, 0.0929912, 0.00102563, 0.000176756, 0.00184695, 0.0261891, 13.4388, 0.395591, 0.0829676, 0.000898367, 0.0907073, 0.000483117, 0.0306672, 0.166496, 0.000362649, 0.399909, 0.00273783, 0.0227258, 3.30817e-6, 8.66783e-6, 9.05333e-6, 2.4177e-5, 0.000184217, 2.80325e-6, 3.56209, 1.07825, 3.95194, 0.447834, 0.0389113, 0.110297, 2.13567, 0.27307, 3.36894, 0.368404, 0.000296116, 0.0415469, 4.95713e-6, 2.47911e-5, 2.06936e-5, 3.39401e-5, 0.0908003, 3.87974e-6, 9.97915, 9.99163, 0.288515, 0.632015, 2.12838, 0.320644, 6.53941, 0.325873]
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
        legendfont = Plots.font("Helvetica", 18),
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

    x = ["BEZ235", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib", "Trametinib"]
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


function plot_dead_acc(concs, drugs, efcs)

    g = zeros(182, 8, 6, 8) # total
    t = LinRange(0.0, 96.0, 182)
    d = zeros(182, 6, 6) # datapoints x concs x drugs
    ls = [1, 4, 5, 6, 7, 8]
    for i = 1:6
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
    intg = zeros(182, 6, 6)
    for i = 1:6
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
            legendfont = Plots.font("Helvetica", 18),
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
        return p1
    end

    p = [single_dead(intg, j) for j=1:6]
    return p

end


function figure70()
    ENV["GKSwstype"]="nul"
    tensor, names, concs, conds, _, _ = DrugResponseModel.hcc_all()
    ps = [9.20945, 2.12023, 2.98291, 0.430111, 0.09537, 0.0956342, 0.242984, 0.0828915, 1.71385, 0.299598, 0.000505937, 0.0647114, 5.01559e-6, 1.49885e-5, 0.00544226, 4.35115e-6, 3.03577e-5, 2.15171e-5, 13.9863, 2.77415, 3.06754, 3.09629, 0.129791, 0.002272, 2.12184, 0.461891, 0.362844, 0.0293403, 0.000990191, 0.00063532, 2.98906e-5, 3.40495e-6, 0.000103298, 5.23445e-7, 7.58079e-6, 0.00361932, 1.11528, 1.93197, 3.82999, 3.8369, 7.05901e-6, 3.52973e-5, 0.714193, 0.0650093, 0.521932, 0.0465728, 0.00111718, 6.65118e-6, 0.00254141, 0.0378996, 9.51571e-5, 7.83845e-7, 0.000242097, 0.00102866, 1.5658, 3.63417, 3.88274, 3.97669, 0.0956175, 0.64807, 2.13597, 0.396049, 1.0315, 0.395203, 0.00162117, 0.00053917, 0.000125216, 0.0929912, 0.00102563, 0.000176756, 0.00184695, 0.0261891, 13.4388, 0.395591, 0.0829676, 0.000898367, 0.0907073, 0.000483117, 0.0306672, 0.166496, 0.000362649, 0.399909, 0.00273783, 0.0227258, 3.30817e-6, 8.66783e-6, 9.05333e-6, 2.4177e-5, 0.000184217, 2.80325e-6, 3.56209, 1.07825, 3.95194, 0.447834, 0.0389113, 0.110297, 2.13567, 0.27307, 3.36894, 0.368404, 0.000296116, 0.0415469, 4.95713e-6, 2.47911e-5, 2.06936e-5, 3.39401e-5, 0.0908003, 3.87974e-6, 9.97915, 9.99163, 0.288515, 0.632015, 2.12838, 0.320644, 6.53941, 0.325873]
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

    # cell deaths
    deathContG1 = zeros(4, 6)
    deathEC50G1 = zeros(4, 6)
    deathContG1[1, :] .= (efcs[9, 1, :]) ./ (efcs[9, 1, :] .+ efcs[1, 1, :])
    deathEC50G1[1, :] .= (ec50[9, :]) ./ (ec50[9, :] .+ ec50[1, :])
    for i = 2:4
        deathContG1[i, :] .= (1 .- deathContG1[i - 1, :]) .* (efcs[i + 8, 1, :]) ./ (efcs[i + 8, 1, :] .+ efcs[i, 1, :])
        deathEC50G1[i, :] .= (1 .- deathEC50G1[i - 1, :]) .* (ec50[i + 8, :]) ./ (ec50[i + 8, :] .+ ec50[i, :])
    end
    deathContG2 = zeros(4, 6)
    deathEC50G2 = zeros(4, 6)
    deathContG2[1, :] .= (efcs[13, 1, :]) ./ (efcs[13, 1, :] .+ efcs[5, 1, :])
    deathEC50G2[1, :] .= (ec50[13, :]) ./ (ec50[13, :] .+ ec50[5, :])
    for i = 14:16
        deathContG2[i - 12, :] = (1 .- deathContG2[i - 13, :]) .* (efcs[i, 1, :]) ./ (efcs[i, 1, :] .+ efcs[i - 8, 1, :])
        deathEC50G2[i - 12, :] = (1 .- deathEC50G2[i - 13, :]) .* (ec50[i, :]) ./ (ec50[i, :] .+ ec50[i - 8, :])
    end


    p0 = Plots.plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p1 = plot_fig1(concs[1], Gshort[:, :, 1, 1], gmshort[:, :, 1, 1], "Dynamical Model Fits - BEZ235", "G1", "A", :PuBu_6, t)
    p2 = plot_fig1(concs[1], Gshort[:, :, 1, 2], gmshort[:, :, 1, 2], "Dynamical Model Fits - BEZ235", "S-G2", "B", :PuBu_6, t)
    p3 = plot_fig1(concs[2], Gshort[:, :, 2, 1], gmshort[:, :, 2, 1], "Dynamical Model Fits - Doxorubicin", "G1", "C", :PuBu_6, t)
    p4 = plot_fig1(concs[2], Gshort[:, :, 2, 2], gmshort[:, :, 2, 2], "Dynamical Model Fits - Doxorubicin", "S-G2", "D", :PuBu_6, t)
    p5 = plot_fig1(concs[3], Gshort[:, :, 3, 1], gmshort[:, :, 3, 1], "Dynamical Model Fits - Gemcitabine", "G1", "E", :PuBu_6, t)
    p6 = plot_fig1(concs[3], Gshort[:, :, 3, 2], gmshort[:, :, 3, 2], "Dynamical Model Fits - Gemcitabine", "S-G2", "F", :PuBu_6, t)
    p7 = plot_fig1(concs[4], Gshort[:, :, 4, 1], gmshort[:, :, 4, 1], "Dynamical Model Fits - Paclitaxel", "G1", "G", :PuBu_6, t)
    p8 = plot_fig1(concs[4], Gshort[:, :, 4, 2], gmshort[:, :, 4, 2], "Dynamical Model Fits - Paclitaxel", "S-G2", "H", :PuBu_6, t)
    p9 = plot_fig1(concs[5], Gshort[:, :, 5, 1], gmshort[:, :, 5, 1], "Dynamical Model Fits - Palbociclib", "G1", "I", :PuBu_6, t)
    p10 = plot_fig1(concs[5], Gshort[:, :, 5, 2], gmshort[:, :, 5, 2], "Dynamical Model Fits - Palbociclib", "S-G2", "J", :PuBu_6, t)
    p11 = plot_fig1(concs[6], Gshort[:, :, 6, 1], gmshort[:, :, 6, 1], "Dynamical Model Fits - Trametinib", "G1", "K", :PuBu_6, t)
    p12 = plot_fig1(concs[6], Gshort[:, :, 6, 2], gmshort[:, :, 6, 2], "Dynamical Model Fits - Trametinib", "S-G2", "L", :PuBu_6, t)
    p = plot_dead_acc(concs, names, efcs)

    figure1 = Plots.plot(p1, p2, p[1], p3, p4, p[2], p5, p6, p[3], p7, p8, p[4], p9, p10, p[5], p11, p12, p[6], size = (1500, 3000), layout = (6, 3))
    Plots.savefig(figure1, "figure70.svg")
end
