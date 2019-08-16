#### Run "Pkg.test("DrugResponseModel") to run the tests "
using Test, DrugResponseModel, Base.Test

println("Starting tests")

tic()
println("Test 1")
@time @test include("testDDE.jl")
println("Test 2")
@time @test include("testHill.jl")
toc()
