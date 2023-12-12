module SVDTreeCompression

using TSVD
using LinearAlgebra: Diagonal
# using Decomposition: split_view


MatrixOrView = Union{ Matrix, SubArray }


mutable struct Node5
    x
    y
    z
    function Node5(y)
        foo = new()
        foo.y = y
        return foo
    end
end

mutable struct Tree_Node
    matrix :: MatrixOrView
    top_left :: Tuple{Int, Int}
    top_right_child :: Tree_Node
    top_left_child :: Tree_Node
    bottom_left_child :: Tree_Node
    bottom_right_child :: Tree_Node
    is_zeros :: Bool
    is_compressed :: Bool
    u :: MatrixOrView
    s :: Vector
    v :: MatrixOrView
    svd_error :: Float64

    function Tree_Node(matrix, top_left)
        foo = new()
        foo.matrix = matrix
        foo.top_left = top_left
        return foo
    end
end

function create_tree(A::MatrixOrView, top_left::Tuple{Int, Int} = (1, 1), r::Int64 = 1, ϵ::Float64 = 1.0)
    node = Tree_Node(A, top_left)
    node.is_zeros = iszero(A)
    if node.is_zeros
        return node
    elseif min(size(A)...) < 16
        node.is_compressed = false
        return node
    end
    u, s, v = tsvd(A, r + 1)
    if s[r + 1] < ϵ
        node.is_compressed = true
        node.u = u[:, 1:r]
        node.s = s[1:r]
        node.v = v[:, 1:r]  # na pseudokodzie jest tu mnozenie Diagonal(node.s) * transpose(node.v)
        return node
    end

    n, _ = div.(size(A), 2)
    node.top_left_child = create_tree(@views A[begin:n, begin:n], (1, 1), r, ϵ)
    node.top_right_child = create_tree(@views A[begin:n, n+1:end], (1, n+1), r, ϵ)
    node.bottom_left_child = create_tree(@views A[n+1:end, begin:n], (n+1, 1), r, ϵ)
    node.bottom_right_child = create_tree(@views A[n+1:end, n+1:end], (n+1, n+1), r, ϵ)
    return node
end

function count_total_tree_error(node)
    if node.is_zeros
        return 0
    elseif node.is_compressed
        return compare_matrixes(node.matrix, node.u*Diagonal(node.s)*transpose(node.v))
    end

    error = 0
    error += count_total_tree_error(node.top_right_child)
    error += count_total_tree_error(node.top_left_child)
    error += count_total_tree_error(node.bottom_left_child)
    error += count_total_tree_error(node.bottom_right_child)
    return error

end

# TODO: Move to Common.jl or sth
function split_view(matrix::MatrixOrView)
    n, _ = div.(size(matrix), 2)
    A11 = @views matrix[begin:n, begin:n]
    A12 = @views matrix[begin:n, n+1:end]
    A21 = @views matrix[n+1:end, begin:n]
    A22 = @views matrix[n+1:end, n+1:end]

    return A11, A12, A21, A22
end

function get_random_nonzero_matrix(size, zeros_percent = 0)
    matrix = rand(size, size)
    println(size ^ 2 * zeros_percent ÷ 100)
    for i in 1:(size ^ 2 * zeros_percent ÷ 100)
        index = rand(1:size, 2)
        while matrix[index[1], index[2]] == 0
            index = rand(1:size, 2)
        end
        matrix[index[1], index[2]] = 0
    end
    return matrix
end

# function create_tree(A, r, ϵ)
#     ...
# end

function compare_matrixes(mat1, mat2)
    return sum((mat1 - mat2) .^ 2)    
end

# xd = 16
# matrix = get_random_nonzero_matrix(xd, 90)
# root = create_tree(matrix, (1, 1), 4, 5.0)


# n, _ = div.(size(matrix), 2)
# matrix2 = @views matrix[begin:n, begin:n]
# println(matrix)

# matrix = [5 5 6 8 3; 7 4 2 3 4; 3 4 9 4 9; 4 7 9 1 5; 8 4 1 5 7]



# println(matrix)
# println(size(u))
# println(size(s))
# println(size(v))265
# println(u * Diagonal(s))


# minsqrdiff = Inf
# mini = 0
# for i in 1:xd-1
#     global minsqrdiff, mini
#     u, s, v = tsvd(matrix, i)
#     svdmatrix = u*Diagonal(s)*transpose(v)
#     # println(round.(svdmatrix, digits = 5))
#     sqrdiff = compare_matrixes(matrix, svdmatrix)
#     if sqrdiff < minsqrdiff
#         minsqrdiff = sqrdiff
#         mini = i
#     end
#     println(i, " ", round.(sqrdiff, digits = 5))
# end

# println("min ", mini, " ", minsqrdiff)
# println("3.2569094358862474 vs ", minsqrdiff * 1024)


end