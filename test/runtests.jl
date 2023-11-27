using MatrixAlgorithms
using LinearAlgebra
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

@testset "IntegerBinet" begin
    for k in 1:256
        A, B = rand(1:100, k, k), rand(1:100, k, k)
        @test Binet(A, B) == AugmentMatrix(A) * AugmentMatrix(B)
    end
end

@testset "FloatBinet" begin
    for k in 1:256
        A, B = Float64.(rand(k, k)), Float64.(rand(k, k))
        @test isapprox(Binet(A, B), AugmentMatrix(A) * AugmentMatrix(B), atol=10^(-8))
    end
end

@testset "IntegerInverse" begin
    for k in 1:256
        A = rand(1:100, k, k)
        @test Inverse(A) * A == I
    end
end

@testset "FloatInverse" begin
    for k in 1:256
        A = Float64.(rand(k, k))
        @test isapprox(Inverse(A) * A, I, atol=10^(-8))
    end
end

@testset "IntegerDecompositionLU" begin
    for k in 1:256
        A = rand(1:100, k, k)
        L, U = DecompositionLU(A)
        @test L * U == A
    end
end

@testset "FloatDecompositionLU" begin
    for k in 1:256
        A = Float64.(rand(k, k))
        L, U = DecompositionLU(A)
        @test isapprox(L * U, A, atol=10^(-8))
    end
end

@testset "IntegerDeterminant" begin
    for k in 1:256
        A = rand(1:100, k, k)
        determinant = Determinant(A)
        @test determinant == det(A)
    end
end

@testset "FloatDeterminant" begin
    for k in 1:256
        A = Float64.(rand(k, k))
        determinant = Determinant(A)
        @test isapprox(determinant, det(A), atol=10^(-8))
    end
end