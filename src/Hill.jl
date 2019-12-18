"""
This file fits Hill function to the parameters.
"""

hill(p, concentration) = p[2] + ((p[3]-p[2])/(1 + ((p[1])/(concentration))^p[4]))

function residHill(hillParams::Array{Float64,1}, concentrations::Array{Float64,1}, g1::Array{Float64,2}, g2::Array{Float64,2}, g1_0::Array{Float64,1}, g2_0::Array{Float64,1}, nG1::Int, nG2::Int)
    """ This functions takes in hill parameters for all the concentrations and calculates
    DDE parameters, passes them to residual function and based off of these, optimizes the model
    and estimates hill parameters. """
    residues = 0.0

    for ii in 1:length(concentrations)
        # [EC50, left, right, steepness]
        alpha =  hill([hillParams[1], hillParams[2], hillParams[3], hillParams[4]], concentrations[ii])
        beta =   hill([hillParams[1], hillParams[5], hillParams[6], hillParams[4]], concentrations[ii])
        gamma1 = hill([hillParams[7], 0.0, hillParams[8], hillParams[9]], concentrations[ii])
        gamma2 = hill([hillParams[7], 0.0, hillParams[10], hillParams[9]], concentrations[ii])

        # collecting all the DDE model parameters
        pp = [alpha, beta, gamma1, gamma2]
        # calculating the residulas for this set of parameters
        residues += cost(pp, g1_0[ii], g2_0[ii], g1[:,ii], g2[:,ii], nG1, nG2)
    end 
    return residues
end

# optimization function for the Hill model
function optimize_hill(guess::Array{Float64,1}, concentrations::Array{Float64,1}, g1::Array{Float64,2}, g2::Array{Float64,2}, g1_0::Array{Float64,1}, g2_0::Array{Float64,1}, lowEC50::Float64, highEC50::Float64, nG1::Int, nG2::Int)
    """ A function to do the optimization given the lower and upper bound of estimation space. """
    # changing the objective function to be compatible with bboptimize
    residue(hillParams) = residHill(hillParams, concentrations, g1, g2, g1_0, g2_0, nG1, nG2)
    # lower bound
    low = [lowEC50, 1e-5, 1e-5, 1e-5, lowEC50, 1e-5, 1e-5, 1e-5, lowEC50, 1e-5, 1e-5, 1e-5]
    # upper bound
    high = [highEC50, 3.0, 3.0, 3.0, highEC50, 3.0, 3.0, 3.0, highEC50, 3.0, 3.0, 3.0]

    res = bboptimize(residue; SearchRange=collect(zip(low, high)), MaxSteps=60000, TraceMode=:silent)
    return best_fitness(res), best_candidate(res)
end

function getODEparams(p::Array{Float64,1}, concentrations::Array{Float64,1})
    """ A function to convert the estimated hill parameters back to ODE parameters. """
    effects = zeros(4, 8)
    for i in 1:8
        effects[1, i] = hill([p[1], p[2], p[3], p[4]], concentrations[i])
        effects[2, i] = hill([p[1], p[5], p[6], p[4]], concentrations[i])
        effects[3, i] = hill([p[7], 0.0,  p[8], p[9]], concentrations[i])
        effects[4, i] = hill([p[7], 0.0, p[10], p[9]], concentrations[i])
    end
    return effects
end