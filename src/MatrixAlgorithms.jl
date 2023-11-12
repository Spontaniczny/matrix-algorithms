using Revise

module MatrixAlgorithms

export Strassen, SplitMatrix

function SplitMatrix(A)
    n, _ = div.(size(A), 2)
    A11 = A[begin:n, begin:n]
    A12 = A[begin:n, n+1:end]
    A21 = A[n+1:end, begin:n]
    A22 = A[n+1:end, n+1:end]

    return A11, A12, A21, A22
end

function Strassen(A, B)
    if size(A) == (1, 1) || size(B) == (1, 1)
        return A * B
    end

    A11, A12, A21, A22 = SplitMatrix(A)
    B11, B12, B21, B22 = SplitMatrix(B)
    
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

end # module MatrixAlgorithms
