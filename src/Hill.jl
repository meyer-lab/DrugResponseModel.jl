"""
This file fits Hill function to the parameters.
"""

hill(p, concentration) = p[2] + ((p[3]-p[2])/(1 + ((p[1])/(concentration))^p[4]))

""" This functions takes in hill parameters for all the concentrations and calculates
DDE parameters, passes them to residual function and based off of these, optimizes the model
and estimates hill parameters. """
function residHill(hillParams::Array{Float64,1}, concentrations::Array{Float64,1}, g1::Array{Float64,2}, g2::Array{Float64,2}, g1_0::Array{Float64,1}, g2_0::Array{Float64,1})
    residues = 0.0

    params = getODEparams(hillParams, concentrations)
    for ii in 1:length(concentrations)
        # collecting all the DDE model parameters
        pp = [params[1,ii], params[2,ii], params[3,ii], params[4,ii], params[5,ii], params[6,ii], params[7,ii]]
        # calculating the residulas for this set of parameters
        residues += cost(pp, g1_0[ii], g2_0[ii], g1[:,ii], g2[:,ii], pp[6], pp[7])
    end 
    return residues
end

""" Hill optimization function. """
function optimize_hill(lowEC50::Float64, highEC50::Float64, conc_l::Array{Float64,1}, g1::Array{Float64,2}, g2::Array{Float64,2}, g1_0::Array{Float64,1}, g2_0::Array{Float64,1})
    lowEC50 = 50.0
    highEC50 = 350.0
    hillCost(hillParams) = residHill(hillParams, conc_l, g1, g2, g1_0, g2_0)
#     hillCost(hillParams) = fullHillCost(hillParams, conc_l, g1, g2, g1_0, g2_0) 
    low =  [lowEC50, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, lowEC50, 1e-5, 1e-5, 1e-5, 0.0, 1, 1]
    high = [highEC50, 3.0, 3.0, 3.0, 3.0, 3.0, highEC50, 3.0, 3.0, 3.0, 1.0, 60, 60]

    results_ode = bboptimize(hillCost; SearchRange=collect(zip(low, high)),
                                           NumDimensions=length(low),
                                           TraceMode=:verbose,
                                           TraceInterval=50,
                                           MaxSteps=4E4);

    return best_fitness(results_ode), best_candidate(results_ode)
end

""" A function to convert the estimated hill parameters back to ODE parameters. """
function getODEparams(p::Array{Float64,1}, concentrations::Array{Float64,1})

    effects = zeros(7, 8)
    for i in 1:8
        # [EC50, left, right, steepness]
        effects[1, i] = hill([p[1], p[2], p[3], p[4]], concentrations[i])
        effects[2, i] = hill([p[1], p[5], p[6], p[4]], concentrations[i])
        effects[3, i] = hill([p[7], 0.0,  p[8], p[9]], concentrations[i])
        effects[4, i] = hill([p[7], 0.0, p[10], p[9]], concentrations[i])
        effects[5, i] = p[11]
        effects[6, i] = p[12]
        effects[7, i] = p[13]
    end
    return effects
end