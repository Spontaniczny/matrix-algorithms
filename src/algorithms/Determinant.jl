function Determinant(A)
    n, _ = size(A)
    if size(A) != (n, n)
        return "Error"
    elseif n == 1
        return A[1, 1]
    elseif n == 2
        return A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    elseif n == 3
        a1 = A[1, 1] * A[2, 2] * A[3, 3] - A[1, 3] * A[2, 2] * A[3, 1]
        a2 = A[2, 1] * A[3, 2] * A[1, 3] - A[2, 3] * A[3, 2] * A[1, 1]
        a3 = A[3, 1] * A[1, 2] * A[2, 3] - A[3, 3] * A[1, 2] * A[2, 1]
        return a1 + a2 + a3 
    end

    L, U = DecompositionLU(A)
    det = 1
    if L[1, 1] == 1
        for i in 1:n
            det *= U[i, i]
        end
        return det
    end

    for i in 1:n
        det *= U[i, i] * L[i, i]
    end

    return det
end