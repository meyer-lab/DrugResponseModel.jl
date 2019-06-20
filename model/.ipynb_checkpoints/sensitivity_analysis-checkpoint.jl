using DiffEqSensitivity, ForwardDiff, DiffResults

# ODE model

function ODEmodel(du, u, p, t)
    # p = [alpha, beta, gamma1, gamma2, initg1, initg2]
    du[1] = -p[1]*u[1] + 2*p[2]*u[2] - p[3]*u[1]
    du[2] = p[1]*u[1] - p[2]*u[2] - p[4]*u[2]
end
    
tspan = (0.0, 95.5)
p = [1.1419,0.8303,0.7368,0.1588]
u0 = [3.0, 2.0]
prob = ODEProblem(ODEmodel, u0, tspan, p)


function test_odemodel(p)
  _prob = remake(prob;u0=convert.(eltype(p),prob.u0),p=p)
  solve(_prob,Vern9(),save_everystep=false)[end]
end


fd_res = ForwardDiff.jacobian(test_odemodel,p)


res = DiffResults.JacobianResult(u0,p) # Build the results object
DiffResults.jacobian!(res,p) # Populate it with the results
val = DiffResults.value(res) # This is the sol[end]
jac = DiffResults.jacobian(res) # This is dsol/dp

