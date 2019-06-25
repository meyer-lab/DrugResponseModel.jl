""" In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle """
import CSV
using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase, Optim, Plots

# Import the data 
data = CSV.read("..//data//Gem.csv")
total = CSV.read("..//data//Gem_pop.csv")

total_old = total[:,4];
G2_old = data[:,4];

estim_init = [3.0, 2.0]
# rescaling the experimental data
total_new = (estim_init[1] + estim_init[2]) * total_old
G2_new = 0.01*total_new.*G2_old;
G1_new = total_new - G2_new;


##---------------------------- Building the function and residuals -------------------------##
function ODEmodel(du, u, p, t)
    # p = [alpha, beta, gamma1, gamma2, initg1, initg2]
    du[1] = -p[1]*u[1] + 2*p[2]*u[2] - p[3]*u[1]
    du[2] = p[1]*u[1] - p[2]*u[2] - p[4]*u[2]
end

function ODEsolve(par)
    t = LinRange(0.0, 95.5, 192)
    tspan = (0.0, 95.5)
    u0 = [3.0, 2.0]
    prob = ODEProblem(ODEmodel, u0, tspan, par)
    solve(prob, Tsit5())
end

function residuals(par)
    t = LinRange(0.0, 95.5, 192)
    #     res= ((sol(t, idxs=2).u).*100 ./(sol(t, idxs=1).u + sol(t, idxs=2).u)) - G2
    res = zeros(2, length(t))
    sol = ODEsolve(par)
    res[1, :] = sol(t, idxs=1).u - G1_new
    res[2, :] = sol(t, idxs=2).u - G2_new
    return res
end

# first guess
p = [1.1419,0.8303,0.7368,0.1588]
# time range
t = LinRange(0.0, 95.5, 192)
# lower bound and higher bound for parameteres
low = [0.0, 0.0, 0.0, 0.0]
high = [10.0, 10.0, 10.0, 10.0]

result_ode = optimize(residuals, p, Dogleg(), lower = low, upper = high, iterations=10^7)

##------------------------ Solve and plot with estimated parameters ----------------------------## 

params = result_ode.minimizer
total_n = 5.0 * total_old
G2_n = 0.01*total_n.*G2_old
G1_n = total_n - G2_n

sol = ODEsolve(params)

plot(t, sol(t, idxs=2).u, label = "G2 est", title = "Gem Trial 4 ODE", legend=:topleft)
plot!(t, G2_n, label = "G2 new")
plot!(t, (sol(t, idxs=2).u + sol(t, idxs=1).u), label = "total est")
plot!(t, total_n, label = "total new")
savefig("gem_4_ode.png")