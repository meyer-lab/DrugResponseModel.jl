#### Run "Pkg.test("DrugResponseModel") to run the tests "
using Test, DrugResponseModel, Base.Test

println("Starting tests")

include("testDDE.jl")
include("testHill.jl")
