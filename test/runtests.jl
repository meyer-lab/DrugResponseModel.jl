using Test, DrugResponseModel, Profile

println("Starting tests")

include("testODE.jl")
include("testDDE.jl")
include("testHill.jl")
