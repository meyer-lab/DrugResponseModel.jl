import CSV
using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase, Optim, LinearAlgebra

data = CSV.read("DOX.csv")
total = CSV.read("DOX_pop.csv")

total_old = total[:,7]
G2_old = data[:,7]

estim_init = [3.0, 3.4] # [init_g1, init_g2]
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
    solve(prob, MethodOfSteps(Tsit5()), lower = [0.0, 0.0])
end

function resid(pp)
    t = LinRange(0.0, 95.5, 192)
    res = zeros(2, length(t))
    sol = solution(pp)
    res[1,:] = sol(t, idxs=1) - G1_new
    res[2,:] = sol(t, idxs=2) - G2_new
    return res
end

p = [1.078,0.869,0.7141,0.2419,20.0988,20.677,2.834,2.719, 3.0]
lowwer = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
upp = [3.0, 3.0, 3.0, 3.0, 30.0, 30.0, 10.0, 10.0, 10.0]
results = optimize(resid, p, Dogleg(), lower = lowwer, upper = upp)

## ------------------ Solve with new params ----------------------- ##

params = results.minimizer
total_new = (params[7] + params[8]) * total_old
G2_new = 0.01*total_new.*G2_old
G1_new = total_new - G2_new

sol = solution(params)
t = LinRange(0.0, 95.5, 192)

using Plots; 
plot(t, sol(t, idxs=2).u, label = "G2 est", title = "Dox Trial 7 Dogleg()")
plot!(t, G2_new, label = "G2 new")
plot!(t, (sol(t, idxs=2).u + sol(t, idxs=1).u), label = "total est")
plot!(t, total_new, label = "total new")
savefig("DOX7_Dogleg.png")






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
