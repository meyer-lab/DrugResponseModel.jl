# Drug Response Model

![Build](https://github.com/meyer-lab/DrugResponseModel.jl/workflows/Build/badge.svg)
![Test](https://github.com/meyer-lab/DrugResponseModel.jl/workflows/Test/badge.svg)


`PhaseChain` is a Julia package for analyzing drug response in population cell number data with respect to cell cycle phase effects. The model is a system of linear first order differential equations where the parameters of the model reflect the quantified effects of drugs on different cell cycle phases; the rate of phase progression and the rate of cell death. With this model and the inferred parameters from fitting, we predicted the effects of drug combinations.

- [Overview](#Overview)
- [Documentation](#Documentation)
- [Systems Requirements](#system-requirements)
- [Installation Guide](#Installation-Guide)
- [Demo Instructions for Use](#Demo)

# Overview

`PhaseChain` is an open-source Julia package that implements a system of ordinary differential equations for cell counts in G1 and S/G2 phases of the cell cycle. The purpose of this model is to quantify the effects of the drugs on different cell cycle phases and provide predictions of drug combination based on these effects. The input data is in the form of cell counts in G1 and S/G2 cell cycle phases. The model is easily extendable to any number of cell cycle phases that we have cell numbers for. The model has been tested on data from treating AU565 breast cancer cell lines with lapatinib, gemcitabine, doxorubicin, palbociclib, and paclitaxel. 

# Documentation
The `docs` folder includes a few tutorials for getting started with the package. All the functions should have a docstring explaining the purpose of the function, as well as the inputs, outputs, and the type of the variables used.

# System Requirements
## Hardware requirements
`PhaseChain` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements

### OS requirements
This package is supported for *macOS* and *Linux*. The package has been tested on the following systems:
- macOS: Mojave (10.14.1)
- Linux: Ubuntu 20.04

### Julia dependencies
The following list shows the Julia packages used in `PhaseChain`.
`BlackBoxOptim`
`CSV`
`Calculus`
`DSP`
`DataFrames`
`DelimitedFiles`
`Distributions`
`GR`
`LinearAlgebra`
`NumericalIntegration`
`Measures`
`Plots`
`Statistics`
`StatsBase`
`StatsPlots`
`XLSX` 

# Installation Guide

### Clone from GitHub
```
git clone https://github.com/meyer-lab/DrugResponseModel.jl.git
```
It takes a few seconds to clone the repository.

# Demo and Instructions for Use

The folder `src` includes all the functions for importing, fitting, and analyzing the data with the model. Functions in `importData.jl` are used for importing the data. The core of the model, jacobian matrix of the ODE system is in the `ODEmodel.jl`. 

```
make output/figure4.svg
```

To run the unit tests:

```
make test
```

To make the manuscript:

```
make output/manuscript.html 
```

#### Creating synthetic data and fitting the model

The following shows how to create a 2-state synthetic lineage of cells with cell fate and cell lifetime observations, fit them to the model, and output the corresponding transition matrix, initial probability matrix, and the estimated parameters for the distribution of each state.

```
import numpy as np
from lineage.states.StateDistributionGamma import StateDistribution
from lineage.LineageTree import LineageTree

# pi: the initial probability vector
pi = np.array([0.6, 0.4], dtype="float")
# This means that the first cell in our lineage in generation 1
# has a 60% change of being state 0 and a 40% chance of being state 1.
# The values of this vector have to add up to 1 because of the
# Law of Total Probability.

# T: transition probability matrix
T = np.array([[0.75, 0.25],
              [0.25, 0.75]], dtype="float")

# State 0 parameters "Resistant"
bern_p0 = 0.99 # a cell fate parameter (Bernoulli distribution)
gamma_a0 = 7 # the shape parameter of the gamma distribution, corresponding to the cell cycle duration
gamma_scale0 = 7 # the scale parameter of the gamma distribution, corresponding to the cell cycle duration

# State 1 parameters "Susceptible"
bern_p1 = 0.88 # a cell fate parameter (Bernoulli distribution)
gamma_a1 = 7 # the shape parameter of the gamma distribution, corresponding to the cell cycle duration
gamma_scale1 = 1 # the scale parameter of the gamma distribution, corresponding to the cell cycle duration

state_obj0 = StateDistribution(bern_p0, gamma_a0, gamma_scale0)
state_obj1 = StateDistribution(bern_p1, gamma_a1, gamma_scale1)

E = [state_obj0, state_obj1]

# creating the synthetic lineage of 15 cells, given the state distributions, transition probability, and the initial probability vector.
lineage = LineageTree.init_from_parameters(pi, T, E, desired_num_cells=2**9 - 1)
```

Now that we have created the lineages as Python objects, we use the following function to fit this data into the model.

```
from lineage.Analyze import Analyze

X = [lineage] # population just contains one lineage
tHMMobj, pred_states_by_lineage, LL = Analyze(X, 2) # find two states

# Estimating the initial probability vector
print(tHMMobj.estimate.pi)

# Estimating the transition probability matrix
print(tHMMobj.estimate.T)

# The total log likelihood
print(LL)

for state in range(lineage.num_states):
    print("State {}:".format(state))
    print("       estimated state:", tHMMobj.estimate.E[state].params)
    print("original parameters given for state:", E[state].params)
    print("\n")
```

#### Importing the experimental data and fitting the model

```
import numpy as np
from lineage.LineageInputOutput import import_exp_data
from lineage.states.StateDistributionGaPhs import StateDistribution
from lineage.LineageTree import LineageTree
from lineage.Analyze import run_Analyze_over

desired_num_states = 2 # does not make a difference what number we choose for importing the data.
E = [StateDistribution() for _ in range(desired_num_states)]

# Importing only one of the replicates of control condition
control1 = [LineageTree(list_of_cells, E) for list_of_cells in import_exp_data(path=r"lineage/data/heiser_data/new_version/AU00601_A5_1_V5.xlsx")]

output = run_Analyze_over([control1], 2, atonce=False)
```
To find the most likely number of states, we can calculate the BIC metrc for 1,2,3,... number of states and find out the likelihoods.
The following calculates the BIC for 2 states, as we chose in the `run_Analyze_over` above.

```
BICs = np.array([oo[0].get_BIC(oo[2], 75, atonce=True)[0] for oo in output])
```

The output of fitting could be the transition matrix:
```
np.array([oo[0].estimate.T for oo in output])
```

initial probability matrix:
```
np.array([oo[0].estimate.pi for oo in output])
```

the assigned cell states lineage by lineage:
```
np.array([oo[1] for oo in output])
```

the distribution parameters for each state:
```
for state in range(2):
    print("State {}:".format(state))
    print("       estimated state:", output[0][0].estimate.E[state].params)
    print("\n")
```

Depending on the number of cells and lineages being used for fitting, the run time for `Analyze` and other similar functions that run the fitting, could takes minutes to hours.

# License
This project is covered under the MIT License.
