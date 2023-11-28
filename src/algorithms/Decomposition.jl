module Decomposition

function DecompositionLU(A)
    if size(A) == (1, 1)
        return hcat(1), hcat(A[1, 1])
    end

    A11, A12, A21, A22 = SplitMatrix(AugmentMatrix(A))
    L11, U11 = DecompositionLU(A11)
    U11_rev = Inverse(U11)
    L21 = Strassen(A21, U11_rev)
    L11_rev = Inverse(L11)
    U12 = Strassen(L11_rev, A12)
    S = A22 - Strassen(L21, U12)
    L22, U22 = DecompositionLU(S)

    L = [L11 zeros(size(L11)); L21 L22]
    U = [U11 U12; zeros(size(U11)) U22]
    return L, U
end

end # Decomposition
