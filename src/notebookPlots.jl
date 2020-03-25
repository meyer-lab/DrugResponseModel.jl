import Pkg;
Pkg.instantiate()
using DrugResponseModel

######################### Lapatinib ##########################
# import data from the path
conc_l, pop, g2, g1 = setup_data("lapatinib")

# costt, ps = optimize_hill(conc_l, g1, g2, maxstep = 4E4)
ps = [68.856, 2.17748, 0.0233308, 0.942434, 1.5169, 2.98134, 0.0181415, 0.0359511, 0.531034, 48.6838, 52.3097, 1.66619, 4.44761]
effects = getODEparams(ps, conc_l)

g0 = g1[1] + g2[1]
T = 170
plotGradient(effects, conc_l, g0, T)

ODEplot_all(effects, g1, g2, pop, conc_l)

ODEplot_allPerc(effects, g1, g2, pop, conc_l)

plot_parameters(conc_l, effects)

result, paramRanges = allSensitivity(ps, conc_l, g1, g2)

using Plots;
pl = [plotUnitSensitivity(paramRanges[:, i], result[:, i], ps[i], i) for i=1:11]
plot(pl...)
plot!(size=(1200, 800))


######################### Doxorubicin ##########################

conc_l, pop, g2, g1 = setup_data("doxorubicin")

# costt, ps = optimize_hill(conc_l, g1, g2, maxstep = 4E4)
ps = [134.546, 0.100021, 0.0385798, 0.814924, 2.93535, 0.229024, 0.0875672, 0.168648, 0.631369, 2.37796, 88.8261, 4.03489, 2.42543]
effects = getODEparams(ps, conc_l)

g0 = g1[1] + g2[1]
T = 170
plotGradient(effects, conc_l, g0, T)

ODEplot_all(effects, g1, g2, pop, conc_l)

ODEplot_allPerc(effects, g1, g2, pop, conc_l)

plot_parameters(conc_l, effects)

result, paramRanges = allSensitivity(ps, conc_l, g1, g2)

pl = [plotUnitSensitivity(paramRanges[:, i], result[:, i], ps[i], i) for i=1:11]
plot(pl...)
plot!(size=(1200, 800))


######################### Gemcitabine ##########################

conc_l, pop, g2, g1 = setup_data("gemcitabine")

# costt, ps = optimize_hill(conc_l, g1, g2, maxstep = 4E4)
ps = [10.8326, 1.00369, 1.03549, 2.19808, 0.444854, 1.17305, 0.00923942, 0.122813, 0.483379, 20.7231, 15.4483, 6.60142, 4.18506]
effects = getODEparams(ps, conc_l)

g0 = g1[1] + g2[1]
T = 170
plotGradient(effects, conc_l, g0, T)

ODEplot_all(effects, g1, g2, pop, conc_l)

ODEplot_allPerc(effects, g1, g2, pop, conc_l)

plot_parameters(conc_l, effects)

result, paramRanges = allSensitivity(ps, conc_l, g1, g2)

using Plots;
pl = [plotUnitSensitivity(paramRanges[:, i], result[:, i], ps[i], i) for i=1:11]
plot(pl...)
plot!(size=(1200, 800))

######################### Paclitaxel ##########################

conc_l, pop, g2, g1 = setup_data("paclitaxel")

# costt, ps = optimize_hill(conc_l, g1, g2, maxstep = 4E4)
ps = [3.62825, 2.43233, 0.521699, 2.58588, 2.34811, 0.818595, 0.0894338, 0.0482111, 0.469062, 47.4889, 69.1083, 3.90522, 1.12117]
effects = getODEparams(ps, conc_l)

g0 = g1[1] + g2[1]
T = 170
plotGradient(effects, conc_l, g0, T)

ODEplot_all(effects, g1, g2, pop, conc_l)

ODEplot_allPerc(effects, g1, g2, pop, conc_l)

plot_parameters(conc_l, effects)

result, paramRanges = allSensitivity(ps, conc_l, g1, g2)

using Plots;
pl = [plotUnitSensitivity(paramRanges[:, i], result[:, i], ps[i], i) for i=1:11]
plot(pl...)
plot!(size=(1200, 800))