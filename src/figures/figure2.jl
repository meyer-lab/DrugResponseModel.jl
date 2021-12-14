""" Plotting the model and data, fitted separately. """

p1 = [4.559212089321431, 1.534630012729537, 6.577953772895624e-6, 0.06015772342568577, 0.08763827716769973, 0.00020696725899275607, 9.900279403817049, 0.3895282713938266, 9.959254675156775, 0.5739737111306142, 0.004501075336460468, 3.1143477217414385e-5, 8.929252102812766e-5, 0.03362498138733098, 0.5007420785015034, 4.599678610545463e-5, 0.0006927823749026561, 0.08033217698520213, 0.056777074388610224, 0.3751064145733935, 0.2815487185584501, 0.3133119169277199, 9.9149790391101, 0.22395916291838106, 9.958324478644093, 9.993999643666973]
p2 = [16.341824946537763, 2.026731958767472, 0.0697586666395713, 0.03769350082536218, 4.732187750376674, 7.917021168287058, 0.3333119356739554, 0.5620294847302345, 9.883686496866622, 9.554838590592949, 0.00338312447921167, 6.147485683149702e-5, 0.0020849651014483335, 0.03532461472377865, 0.004174248975441355, 0.000581256576122533, 0.00026443722232743876, 0.00073376493659534, 9.661831464990229, 0.041129653562686275, 6.869928282102248, 9.483439266169263, 1.4384415187963013, 0.25098501858957456, 9.99559870797076, 9.580143234347688]
p3 = [45.303197545038934, 1.2284309708753, 4.586277547576495, 0.08221921500811005, 0.030923517833211366, 0.3543459386188055, 4.793098625407479, 0.09160401097241552, 4.898741445479735, 0.48658901522898157, 0.00022863841273188982, 3.30717803975358e-5, 1.821899608686722e-5, 4.1323731286719586e-5, 3.1264233954426284e-5, 7.669602385064842e-5, 0.00030337701789132565, 1.2505057782039295e-5, 9.100824309358615, 0.06445751543800647, 0.8656472659383468, 0.27267041534515474, 9.00822037946992, 3.2438791851194653, 5.628164000062103, 0.21050002109104654]
p4 = [465.34625895407356, 2.8337575407378157, 0.08815792513881067, 0.00032883673614319246, 0.004409257403826871, 0.0012518035841512668, 0.33935026914315547, 0.11999892655720527, 1.0501003621090572, 0.5283154636708538, 0.0018143575846060016, 0.01028141316003785, 0.00017878653360945074, 0.0003176802695940445, 0.0003409110851859134, 0.0002491799267869462, 0.009192964268932681, 0.08653656423420368, 0.049391677870604785, 1.5578476759664912, 4.212716165648628, 0.7538437321424133, 7.355358824390203, 2.8421112672082742, 6.619202701159445, 0.2832438493675962]
p5 = [4.900117192623844, 3.907307812428058, 7.241249786823434, 3.1310406893293785, 0.05309849780059287, 0.1380812706380845, 0.05060160472106309, 9.502742066815083, 6.293068224419252, 0.3574471724746114, 0.005338307254828464, 0.0005251038158945597, 4.816902555792851e-5, 2.2497335926667676, 0.002344012168119355, 0.13799970974195316, 0.004414932758367093, 2.9994728472344937e-5, 8.071235872970576, 3.319308815631274, 0.04634269318321133, 8.142580628126892, 0.7123758511538437, 9.906671735027416, 7.018226896985317, 0.36895399085786873]
p6 = [20.424361957329644, 7.598799352261871, 4.155822344229994e-5, 9.445106484982402, 0.027600270736820226, 1.0941197185427223e-5, 6.056972009095306, 7.660130560235781, 0.5369911542697333, 9.998953336035024, 0.06544544289234865, 0.00015174925279038204, 1.2222589873691756e-6, 0.1657441073928077, 0.013828516846045934, 1.123856187132046, 8.911018291135053e-8, 2.64598156158466e-5, 9.999091991412515, 9.75278767765653, 0.03574950674070614, 1.9312162247741684, 9.974662451395005, 9.99351614605866, 0.21954643332238954, 9.999576346521096]
p7 = [0.004973563575347066, 0.10869518657878127, 9.655249346417111, 0.20684808617236153, 0.007690942758531545, 0.2925557959172801, 9.255433337079875, 1.1566129038803066, 0.3274270320278303, 9.827811891574829, 0.00026550017926060506, 0.00400645996901929, 1.1815250363681792e-6, 3.801718349030178e-5, 0.009295088510682309, 0.0004725924310664338, 0.0016113346532294506, 3.12040407123324e-5, 9.958553104492378, 0.04974686104460528, 0.16595230987547438, 0.5369258832356736, 9.289079935937753, 1.1475229621353118, 0.29463634755433965, 9.909359648376803]
p8 = [298.1705678945283, 1.0434564928768324, 0.0002947752846202261, 0.01917020379323064, 7.128820922900845, 0.6833828633983944, 1.3609762742730358, 0.449141637480075, 0.0012542701620185022, 7.388056513336777, 0.026640831446807976, 0.0006790750253830614, 0.5420367389809053, 0.001688341285884588, 0.0004914946237221917, 0.0035418398362248656, 0.09281727520227684, 0.007690856499496632, 7.430311109080693, 0.040863623054482885, 7.322497013931022, 4.690535599281974, 9.58921940940317, 0.22842037607663773, 9.976233553010186, 7.9261447459022785]
p9 = [0.3611827688595246, 0.4397775048024337, 0.11464837672661107, 0.5644655861016754, 0.0026815973004683307, 0.049023679880638854, 3.6761932403582223, 0.40102538761731654, 0.7879239787173435, 9.58646484949174, 3.3341003506554334e-6, 0.00026927978372880357, 0.12763053410783617, 0.0004394629998460365, 9.812050170594746, 9.91606969952789, 4.574197356611854e-5, 0.000555125195411977, 0.04242999241033863, 2.397346001723304, 9.991228585711617, 0.2140129784460572, 4.4543585110030595, 0.39973653093104156, 0.6626399808008144, 9.602894783458693]
p10 = [1468.1056649357215, 2.062089386473211, 0.0021395154032257477, 0.12194842250450827, 0.10991299753964458, 0.6952097118949567, 0.5153934257256693, 0.2620709489284564, 0.010497796821914405, 0.12781635158253038, 0.019812988491970246, 9.414176006326523e-5, 0.001399115209659219, 0.029212515360649505, 5.940578419998891, 2.7327368437275155, 0.008523112333899555, 0.0001420181799657329, 0.5890247413706359, 0.047363083087278136, 0.17634083480006932, 3.393646479002377, 8.863714635272594, 2.2100620138667693, 1.5485281687916483, 0.3664678956591239]
p11 = [3307.2256210598666, 2.825862726889699, 4.261649468429277, 0.19086874320087593, 0.0002088647472333608, 6.151729716215512, 0.6462763902257316, 0.0056944742928132395, 0.31145097527115617, 0.07686550432236511, 0.019940679107192792, 2.106447231956516e-5, 0.003039857098503368, 0.0051940483993526835, 0.09079015881404708, 0.0015977300204902346, 0.0003543178067210582, 0.00953275658601001, 9.79286322569489, 0.07066763207708572, 0.13645503054696045, 9.820276879179776, 8.371531448985223, 0.828331090228319, 0.2809518810892297, 7.8461734618665675]
ps = hcat(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11)
t = LinRange(0.0, 95, 189)

function Eachdrug_sim(G_sim, G_data, condition, g1g2, drug_name)
    df = [DataFrame(x=t, y=G_sim[:, i], y2=G_data[:, i], conc=condition[i]) for i = 1:8]
    DF = vcat(df...)
    p = Gadfly.plot(DF, layer(x="x", y="y", Geom.line, color="conc"), 
                        layer(x="x", y="y2", Geom.line, color="conc", Theme(line_style=[:dash])),
                        Guide.xlabel("time [hr]"), Guide.ylabel("$g1g2 Cell #"), Coord.Cartesian(ymin=-0.05,ymax=2.0), Guide.title(drug_name))
    return p
end

function figure2()
    setGadflyTheme()

    g, c = DrugResponseModel.import_data()
    newg = DrugResponseModel.trim_data(g, c)
    tensor, conditions = DrugResponseModel.form_tensor(newg, c)
    tens = mean(tensor, dims=2)[:, 1, :, :, :]
    pp = []
    for (i, drug) in enumerate(drugs)
        _, _, conc = DrugResponseModel.load_newData(drug)
        pODE = DrugResponseModel.getODEparams(ps[:, i], conc)
        Gs = zeros(189, 8, 2)
        for i=1:8
            Gs[:, i, 1], Gs[:, i, 2], _ = DrugResponseModel.predict(pODE[:, i, 1], pODE[:, 1, 1], t)
        end
        push!(pp, Eachdrug_sim(Gs[:, :, 1], tens[1, :, :, i], conditions[i], "G1", drug))
        push!(pp, Eachdrug_sim(Gs[:, :, 2], tens[2, :, :, i], conditions[i], "S/G2", drug))
    end

    print(length(pp))
    pl = plotGrid((6, 4), [pp..., nothing, nothing];)
    return draw(SVG("figure2.svg", 22inch, 18inch), pl)
end
