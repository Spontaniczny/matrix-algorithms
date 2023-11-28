module Inverse
export inverse, inverse!

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

function inverse(matrix::Matrix)::Matrix
    @assert ispow2(size(matrix, 1))
    matrix |> copy |> inverse!
end

function inverse!(matrix::MatrixOrView)::MatrixOrView
    n = size(matrix, 1)

    if n <= 2
        return inv(matrix)
    end

    A11, A12, A21, A22 = split_view(matrix)
    A11_inverse = inverse!(A11)

    A22 .-= A21 * A11_inverse * A12
    A22_inverse = inverse!(A22)

    A11_inverse_copy = copy(A11_inverse)
    A11_inverse .=  A11_inverse_copy + A11_inverse_copy * A12 * A22_inverse * A21 * A11_inverse_copy
    A12 .= -A11_inverse_copy * A12 * A22_inverse
    A21 .= -A22_inverse * A21 * A11_inverse_copy

    return matrix
end

end # Inverse