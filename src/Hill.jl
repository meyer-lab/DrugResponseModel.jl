"""
This file fits Hill function to the parameters.
"""

hill(p, concentration) = p[2] + ((p[3]-p[2])/(1 + ((p[1])/(concentration))^p[4]))

function residHill(hillParams::Array{Float64,1}, concentrations::Array{Float64,1}, g1::Array{Float64,2}, g2::Array{Float64,2}, g1_0::Array{Float64,1}, g2_0::Array{Float64,1}, nG1::Int, nG2::Int)
    """ This functions takes in hill parameters for all the concentrations and calculates
    DDE parameters, passes them to residual function and based off of these, optimizes the model
    and estimates hill parameters. """
    residues = 0.0

    _, _, params = getODEparams(hillParams, concentrations)
    for ii in 1:length(concentrations)
        # collecting all the DDE model parameters
        pp = [params[1,ii], params[2,ii], params[3,ii], params[4,ii], params[5,ii]]
        # calculating the residulas for this set of parameters
        residues += cost(pp, g1_0[ii], g2_0[ii], g1[:,ii], g2[:,ii], nG1, nG2)
    end 
    return residues
end

""" This helps finding the optimum number of species. """
function fullHillCost(hillParams::Array{Float64,1}, conc_l::Array{Float64,1}, g1::Array{Float64,2}, g2::Array{Float64,2}, g1_0::Array{Float64,1}, g2_0::Array{Float64,1})
    nG1 = floor(hillParams[end-1])
    nG2 = floor(hillParams[end])
    
    return DrugResponseModel.residHill(hillParams[1:end-2], conc_l, g1, g2, g1_0, g2_0, Int(nG1), Int(nG2))
end

""" Hill optimization function. """
function optimize_hill(lowEC50::Float64, highEC50::Float64, conc_l::Array{Float64,1}, g1::Array{Float64,2}, g2::Array{Float64,2}, g1_0::Array{Float64,1}, g2_0::Array{Float64,1})
    lowEC50 = 50.0
    highEC50 = 350.0
    hillCost(hillParams) = fullHillCost(hillParams, conc_l, g1, g2, g1_0, g2_0) 
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

    effects = zeros(5, 8)
    for i in 1:8
        # [EC50, left, right, steepness]
        effects[1, i] = hill([p[1], p[2], p[3], p[4]], concentrations[i])
        effects[2, i] = hill([p[1], p[5], p[6], p[4]], concentrations[i])
        effects[3, i] = hill([p[7], 0.0,  p[8], p[9]], concentrations[i])
        effects[4, i] = hill([p[7], 0.0, p[10], p[9]], concentrations[i])
        effects[5, i] = p[11]
    end
    nG1 = p[end-1]
    nG2 = p[end]
    return Int(floor(nG1)), Int(floor(nG2)), effects
end