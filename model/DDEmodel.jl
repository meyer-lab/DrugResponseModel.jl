import CSV
data = CSV.read("C:\\Users\\farna\\Desktop\\Gem.csv")
G2 = data[:,3]

using Optim, LeastSquaresOptim, DelayDiffEq
using DifferentialEquations, OrdinaryDiffEq

tau1 = 10.0
tau2 = 20.0
# p = [alpha, beta, gamma1, gamma2]
function G1_G2(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-tau1)[1]) + 2*p[2]*(h(p, t-tau2)[2]) - p[3]*(h(p, t-tau1)[1])
    du[2] = p[1]*(h(p, t-tau1)[1]) - p[2]*(h(p, t-tau2)[2]) - p[4]*(h(p, t-tau2)[2])
end

lags = [tau1, tau2]
h(p, t) = ones(2)
t = LinRange(0.0, 95.5, 192)
p = [1.01, 1.0, 0.1, 0.01]
tspan = (0.0, 95.5)
u0 = [1.0, 1.0]

prob = DDEProblem(G1_G2, u0, h, tspan, p; constant_lags = lags)
alg = MethodOfSteps(Tsit5())
# sol = solve(prob, alg)
# using Plots; plot(sol)
cost = build_loss_objective(prob,t,G2,MethodOfSteps(Tsit5()),maxiters=10000)
results = optimize(cost, u0, Dogleg()) # this is the problematic sentence

