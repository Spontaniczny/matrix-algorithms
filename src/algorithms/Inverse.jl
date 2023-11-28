module Inverse
export inverse, inverse!

MatrixOrView = Union{ Matrix, SubArray }

function split_view(matrix::MatrixOrView)::SubArray
    n, _ = div.(size(matrix), 2)
    return @inbounds @views begin 
        matrix[begin:n, begin:n]
        matrix[begin:n, n+1:end]
        matrix[n+1:end, begin:n]
        matrix[n+1:end, n+1:end]
    end
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