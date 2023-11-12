using MatrixAlgorithms
using Test

@testset "IntegerStrassen" begin
    for k in 1:8
        A, B = rand(1:100, 2^k, 2^k), rand(1:100, 2^k, 2^k)
        @test Strassen(A, B) == A * B
    end
end

@testset "FloatStrassen" begin
    for k in 1:8
        A, B = 100*rand(2^k, 2^k), 100*rand(2^k, 2^k)
        @test isapprox(Strassen(A, B), A * B, atol=10^(-5))
    end
end