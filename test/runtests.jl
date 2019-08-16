#### Run "Pkg.test("DrugResponseModel") to run the tests "
using Test, DrugResponseModel

println("Starting tests")

include("testDDE.jl")
include("testHill.jl")
