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
    # Realizing that we needed to add the line below took two people 2 hours to figure out
    # We still have no idea why the inverse method causes side effects
    # We shall dare not ask again
    A11_inv_copy = copy(A11_inv)
    S22_inv = inverse(S22)
    B11 = A11_inv_copy + A11_inv_copy*A12*S22_inv*A21*A11_inv_copy
    B12 = -A11_inv_copy*A12*S22_inv
    B21 = -S22_inv*A21*A11_inv_copy
    B22 = S22_inv
# 
    return [B11 B12; B21 B22]

    return matrix
end

end # Inverse