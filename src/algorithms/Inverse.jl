module Inverse

function Inverse(A)
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

    M1 = A11_rev = Inverse(A11)
    # M2 missing dou to an error 
    M3 = A21_A11_rev = Strassen(A21, M1)
    M4 = Strassen(M3, A12)
    S22 = A22 - M4.res  
    B22 = S22_rev = Inverse(S22)
    M6 = Strassen(A12, S22_rev)
    M7 = Strassen(M6, A21_A11_rev)
    B11 = Strassen(A11_rev, I + M7) 
    A12_rev = Inverse(A12)
    M10 = Strassen(-A11_rev, A12_rev)
    B12 = Strassen(M10, S22_rev)
    B21 = Strassen(-S22_rev, A21_A11_rev)



    return [B11 B12; B21 B22]
end

end # Inverse