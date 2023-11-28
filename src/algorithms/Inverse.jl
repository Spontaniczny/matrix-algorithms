module Inverse
export inverse, inverse!

using LinearAlgebra

# TODO: Move to Common.jl or sth
function split_view(matrix::Matrix)
    n, _ = div.(size(matrix), 2)
    A11 = matrix[begin:n, begin:n]
    A12 = matrix[begin:n, n+1:end]
    A21 = matrix[n+1:end, begin:n]
    A22 = matrix[n+1:end, n+1:end]

    return A11, A12, A21, A22
end

function inverse(matrix::Matrix)::Matrix
    n = size(matrix, 1)

    if n == 2
        return inv(matrix)
    end

    A11, A12, A21, A22 = split_view(matrix)
    A11_inv = inverse(A11)
    S22 = A22 - A21 * A11_inv * A12
    S22_inv = inverse(S22)
    B11 = A11_inv * (I + A12 * S22_inv * A21 * A11_inv)
    B12 = -A11_inv * inverse(A12) * S22_inv
    B21 = - S22_inv * A21 * A11_inv
    B22 = S22_inv

    return [B11 B12; B21 B22]
end

end # Inverse