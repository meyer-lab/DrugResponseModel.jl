import CSV
data = CSV.read("Gem.csv")
total = CSV.read("Gem_pop.csv")

total_old = total[:,7]
G2_old = data[:,7]

estim_init = [4.25, 2.0] # [init_g1, init_g2]
# rescaling the experimental data
total_new = (estim_init[1] + estim_init[2]) * total_old
G2_new = 0.01*total_new.*G2_old
G1_new = total_new - G2_new

using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase, Optim


# This model doesn't assume delay for dying
function G1_G2(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[5])[1]) + 2*p[2]*(h(p, t-p[6])[2]) - p[3]*u[1]
    du[2] = p[1]*(h(p, t-p[5])[1]) - p[2]*(h(p, t-p[6])[2]) - p[4]*u[2] 
end

p = [1.0045, 0.42, 0.01, 0.01, 19.17196, 16.88276, 4.2931205, 1.96889172, 3.0]

function solution(pp)
    lags = [pp[5], pp[6]]
    h(p, t) = pp[9]*ones(2)
    t = LinRange(0.0, 95.5, 192)
    
    tspan = (0.0, 95.5)
    u0 = [pp[7], pp[8]]
    prob = DDEProblem(G1_G2, [pp[7], pp[8]], h, tspan, pp; constant_lags = lags)
    solve(prob, MethodOfSteps(Tsit5()), lower = [0.0, 0.0], upper = [Inf, Inf])
end

function resid(pp)
    res = zeros(2, length(t))
    sol = solution(pp)
    res[1,:] = sol(t, idxs=1).u - G1_new
    res[2,:] = sol(t, idxs=2).u - G2_new
    return res
end

lowwer = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
results = optimize(resid, p, Dogleg(), lower = lowwer)

## ------------------ Solve with new params ----------------------- ##

function ReSolve(results)
    params = results.minimizer
    lags = [params[5], params[6]]
    u0 = [params[7], params[8]]

    tspan = (0.0, 95.5)
    h(p, t) = params[9]*ones(2)
    t = LinRange(0.0, 95.5, 192)
    total_new = (params[7] + params[8]) * total_old
    G2_new = 0.01*total_new.*G2_old
    G1_new = total_new - G2_new

    estim_prob = DDEProblem(G1_G2, u0, h, tspan, params; constant_lags = lags)
    sol = solve(estim_prob, MethodOfSteps(Tsit5()), lower = [0.0, 0.0], upper = [Inf, Inf])
end

sol = ReSolve(results)

using Plots; 
plot(t, sol(t, idxs=2).u, label = "est G2", title = "Trial 7")
plot!(t, (sol(t, idxs=2).u + sol(t, idxs=1).u), label = "est total")
plot!(t, sol(t, idxs=1).u, label = "G1 est")
savefig("7.png")
