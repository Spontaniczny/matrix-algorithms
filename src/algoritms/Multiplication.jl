module Multiplication

using PaddedViews
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
end # Multiplication