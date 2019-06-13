import CSV
using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase

data = CSV.read("Gem.csv")
total = CSV.read("Gem_pop.csv")

# we could choose the columns of "total" and "data" from [2:8]
total_old = total[:,5]
G2_old = data[:,5]

estim_init = [4.25, 2.0] # [init_g1, init_g2] 
# rescaling the experimental data
total_new = (estim_init[1] + estim_init[2]) * total_old
G2_new = 0.01*total_new.*G2_old
G1_new = total_new - G2_new


# This model doesn't assume delay for dying
function G1_G2(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[5])[1]) + 2*p[2]*(h(p, t-p[6])[2]) - p[3]*u[1]
    du[2] = p[1]*(h(p, t-p[5])[1]) - p[2]*(h(p, t-p[6])[2]) - p[4]*u[2] 
end

p = [1.0045,-0.418798,7.47465,0.3212957,19.17196,16.88276,4.2931205,1.96889172]

function resid(pp)
    lags = [pp[5], pp[6]]
    h(p, t) = 2*ones(2)
    t = LinRange(0.0, 95.5, 192)
    
    tspan = (0.0, 95.5)
    u0 = [pp[7], pp[8]]
    prob = DDEProblem(G1_G2, [pp[7], pp[8]], h, tspan, pp; constant_lags = lags)

    sol = solve(prob, MethodOfSteps(Tsit5()))
    
    res = zeros(2, length(t))
    res[1,:] = sol(t, idxs=1).u - G1_new
    res[2,:] = sol(t, idxs=2).u - G2_new
    return res
end

results = optimize(resid, p, Dogleg())


## ---------------------- solve given the estimated parameters ----------------------------##

params = results.minimizer

# updating the parameters to solve DDE
lags = [params[5], params[6]]
u0 = [params[7], params[8]]

total_new = (p[7] + p[8]) * total_old
G2_new = 0.01*total_new.*G2_old
G1_new = total_new - G2_new

tspan = (0.0, 95.5)
h(p, t) = 2*ones(2)
t = LinRange(0.0, 95.5, 192)

# solve
estim_prob = DDEProblem(G1_G2, u0, h, tspan, params; constant_lags = lags)
sol = solve(estim_prob, MethodOfSteps(Tsit5()))
total_cells = sol(t, idxs=1).u + sol(t, idxs=2).u

# plot
using Plots; 

plot(t, sol(t, idxs=2).u, label = "est G2")
plot!(t, G2_new, label = "G2")
plot!(t, total_cells, label = "est total")
plot!(t, total_new, label = "total")