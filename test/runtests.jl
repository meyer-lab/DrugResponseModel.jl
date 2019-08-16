#### Run "Pkg.test("DrugResponseModel") to run the tests "
using Test, DrugResponseModel

println("Starting tests")


println("Test 1")
@elapsed @test include("testDDE.jl")
println("Test 2")
@elapsed @test include("testHill.jl")

