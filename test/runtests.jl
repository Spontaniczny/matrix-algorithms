using MatrixAlgorithms
using Test

@testset "IntegerStrassen" begin
    for k in 1:256
        A, B = rand(1:100, k, k), rand(1:100, k, k)
        @test Strassen(A, B) == AugmentMatrix(A) * AugmentMatrix(B)
    end
end

@testset "FloatStrassen" begin
    for k in 1:256
        A, B = Float64.(rand(k, k)), Float64.(rand(k, k))
        @test isapprox(Strassen(A, B), AugmentMatrix(A) * AugmentMatrix(B), atol=10^(-8))
    end
end