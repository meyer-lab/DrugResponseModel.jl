using Test, DrugResponseModel, Profile

println("Start tests")

include("testODE.jl")
include("testDDE.jl")
include("testHill.jl")
