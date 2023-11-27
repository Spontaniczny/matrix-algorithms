
module MatrixAlgorithms

using Revise
using PaddedViews
using LinearAlgebra

export Strassen, SplitMatrix, Strassen, CountedStrassen, Binet, CountedBinet, AugmentMatrix

function SplitMatrix(A)
    n, _ = div.(size(A), 2)
    A11 = A[begin:n, begin:n]
    A12 = A[begin:n, n+1:end]
    A21 = A[n+1:end, begin:n]
    A22 = A[n+1:end, n+1:end]

    return A11, A12, A21, A22
end

function AugmentMatrix(A)
    n, _ = size(A)
    if ispow2(n)
        return A
    end

    return collect(PaddedView(0, A, (nextpow(2, n), nextpow(2, n)), (1, 1)))
end

function CountedStrassen(A, B)
    if size(A) == (1, 1) || size(B) == (1, 1)
        return (res=A * B, add=0, mul=1)
    end

    n, _ = size(A)

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))
    B11, B12, B21, B22 = SplitMatrix(AugmentMatrix(B))

    M1 = CountedStrassen(A11 + A22, B11 + B22)
    M2 = CountedStrassen(A21 + A22, B11)
    M3 = CountedStrassen(A11, B12 - B22)
    M4 = CountedStrassen(A22, B21 - B11)
    M5 = CountedStrassen(A11 + A12, B22)
    M6 = CountedStrassen(A21 - A11, B11 + B12)
    M7 = CountedStrassen(A12 - A22, B21 + B22)

    C11 = M1.res + M4.res - M5.res + M7.res
    C12 = M3.res + M5.res
    C21 = M2.res + M4.res
    C22 = M1.res - M2.res + M3.res + M6.res

    Ms = [M1, M2, M3, M4, M5, M6, M7]
    matrix_addition_cost = div(n, 2)^2
    total_additions = sum(map(M -> M.add, Ms)) + matrix_addition_cost * 18
    total_multiplications = sum(map(M -> M.mul, Ms))

    return (res=[C11 C12; C21 C22],
            add=total_additions,
            mul=total_multiplications)
end

function Strassen(A, B)
    if size(A) == (1, 1) || size(B) == (1, 1)
        return A * B
    end

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))
    B11, B12, B21, B22 = SplitMatrix(AugmentMatrix(B))
    
    M1 = Strassen(A11 + A22, (B11 + B22))
    M2 = Strassen(A21 + A22, B11)
    M3 = Strassen(A11, B12 - B22)
    M4 = Strassen(A22, B21 - B11)
    M5 = Strassen(A11 + A12, B22)
    M6 = Strassen(A21 - A11, B11 + B12)
    M7 = Strassen(A12 - A22, B21 + B22)

    C11 = M1 + M4 - M5 + M7
    C12 = M3 + M5
    C21 = M2 + M4
    C22 = M1 - M2 + M3 + M6

    return [C11 C12; C21 C22]
end

function CountedBinet(A, B)
    if size(A) == (1, 1) || size(B) == (1, 1)
        return (res=A * B, add=0, mul=1)
    end

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))
    B11, B12, B21, B22 = SplitMatrix(AugmentMatrix(B))

    M1 = CountedBinet(A11, B11)
    M2 = CountedBinet(A12, B21)
    M3 = CountedBinet(A11, B12)
    M4 = CountedBinet(A12, B22)
    M5 = CountedBinet(A21, B11)
    M6 = CountedBinet(A22, B21)
    M7 = CountedBinet(A21, B12)
    M8 = CountedBinet(A22, B22)

    C11 = M1.res + M2.res
    C12 = M3.res + M4.res
    C21 = M5.res + M6.res
    C22 = M7.res + M8.res

    n, _ = size(A)
    Ms = [M1, M2, M3, M4, M5, M6, M7, M8]
    matrix_addition_cost = div(n, 2)^2
    total_additions = sum(map(M -> M.add, Ms)) + matrix_addition_cost * 4
    total_multiplications = sum(map(M -> M.mul, Ms))

    return (res=[C11 C12; C21 C22],
            add=total_additions,
            mul=total_multiplications)
end

function Binet(A, B)
    if size(A) == (1, 1) || size(B) == (1, 1)
        return A * B
    end

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))
    B11, B12, B21, B22 = SplitMatrix(AugmentMatrix(B))

    M1 = Binet(A11, B11)
    M2 = Binet(A12, B21)
    M3 = Binet(A11, B12)
    M4 = Binet(A12, B22)
    M5 = Binet(A21, B11)
    M6 = Binet(A22, B21)
    M7 = Binet(A21, B12)
    M8 = Binet(A22, B22)

    C11 = M1 + M2
    C12 = M3 + M4
    C21 = M5 + M6
    C22 = M7 + M8

    return [C11 C12; C21 C22]
end

function CountedReverse(A)

    if size(A) == (1, 1)
        return (res=hcat(1 / A[1, 1]), add=0, mul=1)
    end

    if size(A) == (2, 2)
        a, c, b, d = A
        x = a * d - b * c
        mat = [d/x -b/x; -c/x a/x]
        return (res=mat, add=1, mul=6)  # do we count -b as (-1) * b?
    end

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))

    M1 = A11_rev = CountedReverse(A11)
    # M2 missing dou to an error 
    M3 = A21_A11_rev = CountedStrassen(A21, M1.res)
    M4 = CountedStrassen(M3.res, A12)
    S22 = A22 - M4.res  # add to total_additions
    M5 = B22 = S22_rev = CountedReverse(S22)
    M6 = CountedStrassen(A12, S22_rev.res)
    M7 = CountedStrassen(M6.res, A21_A11_rev.res)
    M8 = B11 = CountedStrassen(A11_rev.res, I + M7.res) # add to total_additions
    M9 = A12_rev = CountedReverse(A12)
    M10 = CountedStrassen(-A11_rev.res, A12_rev.res)
    M11 = B12 = CountedStrassen(M10.res, S22_rev.res)
    M12 = B21 = CountedStrassen(-S22_rev.res, A21_A11_rev.res)


    Ms = [M1, M3, M4, M5, M6, M7, M8, M9, M10, M11, M12]
    total_additions = sum(map(M -> M.add, Ms)) + length(A22) + length(A11)
    total_multiplications = sum(map(M -> M.mul, Ms))

    return (res=[B11 B12; B21 B22],
            add=total_additions,
            mul=total_multiplications)
end 

function Reverse(A)
    if size(A) == (1, 1)
        return hcat(1 / A[1, 1])
    end

    if size(A) == (2, 2)
        a, c, b, d = A
        x = a * d - b * c
        mat = [d/x -b/x; -c/x a/x]
        return mat
    end

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))

    M1 = A11_rev = Reverse(A11)
    # M2 missing dou to an error 
    M3 = A21_A11_rev = Strassen(A21, M1)
    M4 = Strassen(M3, A12)
    S22 = A22 - M4.res  
    B22 = S22_rev = Reverse(S22)
    M6 = Strassen(A12, S22_rev)
    M7 = Strassen(M6, A21_A11_rev)
    B11 = Strassen(A11_rev, I + M7) 
    A12_rev = Reverse(A12)
    M10 = Strassen(-A11_rev, A12_rev)
    B12 = Strassen(M10, S22_rev)
    B21 = Strassen(-S22_rev, A21_A11_rev)



    return [B11 B12; B21 B22]
end

function CounterDecompositionLU(A)
    if size(A) == (1, 1)
        return (res = (hcat(1), hcat(A[1, 1])), add=0, mul=1)
    end

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))
    M1 = CountedDecompositionLU(A11)
    L11, U11 = M1.res
    M2 = U11_rev = CountedReverse(U11)
    M3 = L21 = CountedStrassen(A21, U11_rev.res)
    M4 = L11_rev = CountedReverse(L11)
    M5 = U12 = CountedStrassen(L11_rev.res, A12)
    M6 = CountedStrassen(L21.res, U12.res)
    S = A22 - M6.res  # add to total_additions
    M7 = CountedDecompositionLU(S)
    L22, U22 = M7.res

    Ms = [M1, M2, M3, M4, M5, M6, M7]
    total_additions = sum(map(M -> M.add, Ms)) + length(A22)
    total_multiplications = sum(map(M -> M.mul, Ms))

    L = [L11 zeros(size(L11)); L21 L22]
    U = [U11 U12; zeros(size(U11)) U22]

    return (res=(L, U),
            add=total_additions,
            mul=total_multiplications)
end

function DecompositionLU(A)
    if size(A) == (1, 1)
        return hcat(1), hcat(A[1, 1])
    end

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))
    L11, U11 = DecompositionLU(A11)
    U11_rev = Reverse(U11)
    L21 = Strassen(A21, U11_rev)
    L11_rev = Reverse(L11)
    U12 = Strassen(L11_rev, A12)
    S = A22 - Strassen(L21, U12)
    L22, U22 = DecompositionLU(S)

    L = [L11 zeros(size(L11)); L21 L22]
    U = [U11 U12; zeros(size(U11)) U22]
    return L, U


end
end # module MatrixAlgorithms
