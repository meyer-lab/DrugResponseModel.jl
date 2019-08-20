# Drug Response Model

[![Build Status](https://travis-ci.com/meyer-lab/DrugResponseModel.jl.svg?branch=master)](https://travis-ci.com/meyer-lab/DrugResponseModel.jl)

This is to show oscillations in the number of cells in G1 and in G2 phase of cell cycle in some cancer cell lines, in response to different concentrations of a number of chemotherapy drugs, including gemcitabine, doxorubicin, paclitaxel, and lapatinib. 


<!-- ### Here is a way to use the DDEmodel

```ruby
include("importData.jl")
include("DDEmodel.jl")
include("plot.jl")

# import data from the path
pop, g2, g1, g1_0, g2_0 = get_data("..//data//lap.csv", "..//data//lap_pop.csv"); # in which:
# pop: population data
# g1, g2: g1 and g2 data
# initial: initial number of cells in g1 and in g2 at time 0

# This is to load the estimated parameters to be used as "initial guess"
param_lap_dde = CSV.read(".//ddeParams//Lapatinib//params_lap_DDE.csv")

# initial guesses for the parameters
lap = convert(Matrix, param_lap_dde[1:7,2:end]); 


# i is the number of the column we are using from the data (# of trial)
i = 6

# initial guess
p  = [0.02798, 0.025502, 21.3481, 10.2881, 0.0001, 0.0001]

# setting lowest delay for tau1 to be half an hour and for tau2 to be 3 hours.
low = [0.015, 0.003, 0.5, 3.0, 7.0, 0.0001, 0.0001]
upp = [0.075, 0.075, 30.0, 30.0, 100.0, 0.05, 0.05]

# Estimating the parameters for trial i
params = optimIt(p, low, upp, i, g1, g2)

# Plotting the long-term prediction along with the data for trial i
plotIt(params, g1, g2, g1_0, g2_0, pop, i, "Lapatinib")

```
<!-- ![Trial 6 for Lapatinib](.png) --> -->



