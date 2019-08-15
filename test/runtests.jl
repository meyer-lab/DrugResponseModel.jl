#### Run "Pkg.test("DrugResponseModel") to run the tests "
using Test, DrugResponseModel

println("Starting tests")

@testset "testing DDE model" begin include("testDDE.jl") end
@testset "testing Hill model" begin include("testHill.jl") end
