import CSV
using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase, Optim, Plots, Statistics, DataFrames

# Import data
data = CSV.read("..//data//lap.csv")
total = CSV.read("..//data//lap_pop.csv")

# delete the extra index column which strats at 0 
deletecols!(data, :Column1)
deletecols!(total, :Column1)

# get the number of data points
len = length(data[:,1])
time = data[:,1]

# just to make sure all have the same size
total_old = total[:,9];
G2_old = data[:,9];

# rescaling the experimental data assuming we have 20 initial cells for each trial
init_cells = 20.0
total_new = init_cells * total_old
G2_new = 0.01*total_new.*G2_old;
G1_new = total_new .- G2_new;
g2_0 = init_cells*(G2_old[1]/100.0)
g1_0 = init_cells*(1- (G2_old[1]/100.0))

plot(G1_new, label = "G1")
plot!(G2_new, label = "G2")
plot!(total_new, label = "total")
print(g1_0, "\n", g2_0)

# param_holder_lap_dde = zeros(7,9);
# param_holder_lap_ode = zeros(4,9);

## --------------------- Define function and residuals ----------------------- ##
# This model doesn't assume delay for dying
function DDEmodel(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[3])[1]) + 2*p[2]*(h(p, t-p[4])[2]) - p[6]*u[1]
    du[2] = p[1]*(h(p, t-p[3])[1]) - p[2]*(h(p, t-p[4])[2]) - p[7]*u[2]
end

function DDEsolve(pp)
    lags = [pp[3], pp[4]]
    h(p, t) = pp[5]*ones(2)
    t = LinRange(0.0, 95.5, 192)

    tspan = (0.0, 95.5)
    u0 = [g1_0, g2_0]
    prob = DDEProblem(DDEmodel, u0, h, tspan, pp; constant_lags = lags)
    solve(prob)
end

function resid(pp)
    t = LinRange(0.0, 95.5, 192)
    res = zeros(2, 192)
    sol = DDEsolve(pp)
    res[1, :] = sol(t, idxs=1).u - G1_new
    res[2, :] = sol(t, idxs=2).u - G2_new
    return res
end

# We could just find the best guess for the first trial and then just substitude the previous estimated parameter 
# to find the next trial as our initial guess (mostl of the time works well!)
p  = [0.00959424, 0.0256688, 20.9196, 18.9652, 8.06589, 0.00227095, 0.0103966]

# setting lowest delay for tau1 to be half an hour and for tau2 to be 3 hours.
low = [0.0001, 0.0001, 0.45, 2.0, 0.0001, 0.0001, 0.0001]
# results_dde = optimize(resid, p, Dogleg(), lower = low)
results_dde = optimize(resid, p, LevenbergMarquardt(), lower = low)

## -------------------- Solve with new params and plot ---------------------- ##

params = results_dde.minimizer
print(params)
# param_holder_lap_dde = CSV.read("params_dox_DDE.csv")
# param_holder_dox_dde[1:7, 4] = params 
# CSV.write("params_dox_DDE.csv",  DataFrame(param_holder_dox_dde))
# [0.0657397, 0.0808217, 0.451049, 36.2012, 9.34548, 0.016664, 0.0193367]

# param_holder_dox_dde[:, 1] = ["alpha", "beta", "tau1", "tau2", "history", "gamma1", "gamma2"]
# CSV.write("params_lap_DDE.csv",  DataFrame(param_holder_lap_dde))
# a function to find the percent of G2/G1


# params = [0.0484957, 0.0297392, 0.92571, 15.1207, 9.84743, 0.00174232, 0.000253351]
# sol = DDEsolve(params)
t = LinRange(0.0, 95.5, 192)

# plot(t, sol(t, idxs=2).u, label = "G1 est", legend =:topleft)
# plot!(t, G1_new, label = "G1")
# plot!(t, sol(t, idxs=1).u, title = "DOX Trial 2 DDE", label = "G2 est")
# plot!(t, G2_new, label = "G2")
# plot!(t, (sol(t, idxs=2).u + sol(t, idxs=1).u), label = "total est")
# plot!(t, total_new, label = "total")
# savefig("DOX_2_DDE.png")

lags = [params[3], params[4]]
t_new = LinRange(0.0, 195.5, 292)
h(p, t_new) = params[5]*ones(2)
tspan_new = (0.0, 195.5)
u0_new = [g1_0, g2_0]
prob_new = DDEProblem(DDEmodel, u0_new, h, tspan_new, params; constant_lags = lags)
solution = solve(prob_new, MethodOfSteps(Tsit5()))

plot(t_new, solution(t_new, idxs=1).u, label = "G1 est", dpi = 150, xlabel = "time [hours]", ylabel = "# of cells")
plot!(t, G1_new, label = "G1", dpi = 150)
plot!(t_new, solution(t_new, idxs=2).u, label = "G2 est", title = "Lap. Trial 9 DDE", legend=:topright, dpi = 150)
plot!(t, G2_new, label = "G2", dpi = 150)
plot!(t_new, (solution(t_new, idxs=2).u + solution(t_new, idxs=1).u), label = "total est", dpi = 150)
plot!(t, total_new, label = "total", dpi = 150)
# savefig("lap_9_dde_long.png")

function g2Tog1(sol)
    t = LinRange(0.0, 95.5, 192)
    return (sol(t, idxs=2).u ./ sol(t, idxs=1).u) .*100
end



# Obsolete!
## ====================== Use Simulated Annealing =========================##

import CSV
using LeastSquaresOptim, DifferentialEquations, DelayDiffEq
using DiffEqBase, Optim, LinearAlgebra, Plots
data = CSV.read("DOX.csv")
total = CSV.read("DOX_pop.csv")

total_old = total[:,7]
G2_old = data[:,7]

estim_init = [1.88, 3.25] # [init_g1, init_g2]
# rescaling the experimental data
total_new = (estim_init[1] + estim_init[2]) * total_old
G2_new = 0.01*total_new.*G2_old
G1_new = total_new - G2_new
 
# This model doesn't assume delay for dying
function G1_G2(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[5])[1]) + 2*p[2]*(h(p, t-p[6])[2]) - p[3]*u[1]
    du[2] = p[1]*(h(p, t-p[5])[1]) - p[2]*(h(p, t-p[6])[2]) - p[4]*u[2] 
end

function solution(pp)
    lags = [pp[5], pp[6]]
    h(p, t) = pp[9]*ones(2)
    t = LinRange(0.0, 95.5, 192)
    tspan = (0.0, 95.5)
    u0 = [pp[7], pp[8]]

    prob = DDEProblem(G1_G2, [pp[7], pp[8]], h, tspan, pp; constant_lags = lags)
    sol = solve(prob, MethodOfSteps(Tsit5()), lower = [0.0, 0.0], upper = [20.0, 20.0])
end

function resid(pp)
    t = LinRange(0.0, 95.5, 192)
    sol = solution(pp)
    res = sol(t, idxs=2) - G2_new
    
#     res = zeros(2, length(t))
#     res[1,:] = sol(t, idxs=1) - G1_new
#     res[2,:] = sol(t, idxs=2) - G2_new
    return norm(res)
end

p = [0.470,0.455,0.490,0.021,17.613,17.626,2.978,2.773,2.923]
lowwer = [0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 0.0, 0.0, 0.0]
upperr = [100.0, 100.0, 100.0, 80.0, 40.0, 40.0, 20.0, 20.0, 100.0]

res = Optim.optimize(resid, lowwer, upperr, p, SAMIN(rt = 0.7), Optim.Options(iterations = 10^8))


# ---------------- Solve with new parameters from samin() ----------------------#

total_nw = (res.minimizer[7] + res.minimizer[8]) * total_old
G2_new = 0.01*total_nw.*G2_old
G1_new = total_nw - G2_new
t = LinRange(0.0, 95.5, 192)

final = solution(res.minimizer)

plot(t, final(t, idxs=2).u, label = "G2 est", title = "Dox Trial 7 samin(rt=0.7)")
plot!(t, G2_new, label = "G2")
plot!(t, (final(t, idxs=2).u + final(t, idxs=1).u), label = "total est")
plot!(t, total_nw, label = "total")
# plot!(t, final(t, idxs=1).u, label = "est G1")
savefig("DOX_one_7_samin_rt7.png")
