module Decomposition

export lu_decompose

using LinearAlgebra: inv

MatrixOrView = Union{ Matrix, SubArray }

# TODO: Move to Common.jl or sth
function split_view(matrix::MatrixOrView)
    n, _ = div.(size(matrix), 2)
    A11 = @views matrix[begin:n, begin:n]
    A12 = @views matrix[begin:n, n+1:end]
    A21 = @views matrix[n+1:end, begin:n]
    A22 = @views matrix[n+1:end, n+1:end]

    return A11, A12, A21, A22
end

function simple_decompose(matrix::MatrixOrView)::Tuple{Matrix, Matrix}
    @assert size(matrix, 1) == 2
    a = matrix[begin, begin]
    b = matrix[begin, end]
    c = matrix[end, begin]
    d = matrix[end, end]
    L = [1 0; (c / a) 1]
    U = [a b; 0 d - b * c / a]
    return L, U
end

# Yes, this implementation is not up to my standards, but time is short 
function lu_decompose(matrix::MatrixOrView)::Tuple{Matrix, Matrix}
    @assert ispow2(size(matrix, 1))
    if size(matrix, 1) == 2 
        return simple_decompose(matrix)
    end

    A11, A12, A21, A22 = split_view(matrix)
    L11, U11 = lu_decompose(A11)
    U11_inverse = inv(U11)
    L21 = A21 * U11_inverse
    L11_inverse  = inv(L11)
    U12 = L11_inverse * A12
    S = A22 - A21 * U11_inverse * U12
    Ls, Us = lu_decompose(S)
    U22 = Us
    L22 = Ls
    L = [L11 zeros(size(L11)); L21 L22]
    U = [U11 U12; zeros(size(L11)) U22]
    return L, U
end

end # Decomposition
