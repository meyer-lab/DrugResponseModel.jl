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
We observed oscillation over time in the phase-specific cell numbers, so we looked into the single cell data. Because of the distribution of cell cycle phase lengths to use "linear chain trick" and created the mean field system of ODEs from the stochastic state transition process of cell division and phase transition. Please refer to the method section in the manuscript for more details.

# Documentation
The `docs` folder includes a tutorial for getting started with the package. All the functions should have a docstring explaining the purpose of the function, as well as the inputs, outputs, and the type of the variables used.

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

The folder `src` includes all the functions for importing, fitting, and analyzing the data with the model. Functions in `importData.jl` are used for importing the data. The core of the model, jacobian matrix of the ODE system is in the `ODEmodel.jl`. Most plotting functions are in the `plot.jl` and optimization for only one data of drug treatment is in the `Hill.jl`, and for fitting all drugs at once, another fitting function considering some simplifying assumptions are located in the `allDrugs.jl`.
The figures created for the manuscript are in the `src/figures` folder. The unittests are in the `test`.


#### Importing data and fitting

The following shows how to import the data and fit into the model in terminal while inside the repository main folder, assuming you have Julia installed.

```
import Pkg; Pkg.instantiate()
Pkg.activate(".")
using DrugResponseModel

# import one replicate of the data for all 5 drug treatments and the concentrations.
concs, _, g1s1, g2s1 = DrugResponseModel.load(189, 1);
```
`concs` is a [8 x 5] matrix, containing the 8 concentrations (including control) for the 5 drugs.
The indexes corresponding to each drug are: 1. lapatinib, 2. doxorubicin, 3. gemcitabine, 4. paclitaxel, and 5. palbociclib.
`g1s1` and `g2s1` are [189, 8, 5] matrices corresponding to G1 and S/G2 cell cycle phases, where they include cell numbers for 189 time points (96 hours) for the 8 concentrations of the 5 drugs.
The dose-response behavior of all drugs are assumed to be Hill-shaped. cell death rates are assumed strictly increasing and the cell cycle progression rate is assumed to be strictly decreasing over the concentrations.

```
# finding the parameters for lapatinib treatment first replicate
cost1, lapatinib_parameters = DrugResponseModel.optimize_hill(concs[:, 1], g1s1[:, :, 1], g2s1[:, :, 1])

# finding the parameters for gemcitabine treatment first replicate
cost2, gemcitabine_parameters = DrugResponseModel.optimize_hill(concs[:, 3], g1s1[:, :, 3], g2s1[:, :, 3])
```

The parameters that are estimated here, are not directly inputted to the ODE system. They are Hill function parameters for each progression and death rate. The following converts these estimated parameters to the ones that correspond to the ODE system.

```
param_lpt = DrugResponseModel.getODEparams(lapatinib_parameters, concs[:, 1])
param_gem = DrugResponseModel.getODEparams(gemcitabine_parameters, concs[:, 3])
```


# License
This project is covered under the MIT License.
