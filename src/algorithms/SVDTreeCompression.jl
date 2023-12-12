module SVDTreeCompression

using TSVD
using LinearAlgebra: Diagonal


MatrixOrView = Union{ Matrix, SubArray }

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
    if node.is_zeros || min(size(node.matrix)...) < 4
    # if node.is_zeros
    #     return node
    # elseif min(size(node.matrix)...) < 16
        return node
    end

    is_ok = false
    counter = 0
    while !is_ok && counter < 100
        try
            u, s, v = tsvd(A, r + 1)
            is_ok = true

        catch e
            x, _ = size(A)
            index = rand(1:x, 2)
            A[index[1], index[2]] += 1e-15
            counter += 1
            
        end
    end
    u, s, v = tsvd(A, r + 1)
    
    if s[r + 1] < ϵ
        node.is_compressed = true
        node.u = u[:, 1:r]
        node.s = s[1:r]
        node.v = v[:, 1:r]  # na pseudokodzie jest tu mnozenie Diagonal(node.s) * transpose(node.v)
        return node
    end
    A11, A12, A21, A22 = split_view(A)
    n, _ = div.(size(A), 2)

    node.top_left_child = create_tree(A11, (1, 1), r, ϵ) # ask @integraledelebesgue why @view does not work here
    node.top_right_child = create_tree(A12, (1, n+1), r, ϵ) # ask @integraledelebesgue why @view does not work here
    node.bottom_left_child = create_tree(A21, (n+1, 1), r, ϵ) # ask @integraledelebesgue why @view does not work here
    node.bottom_right_child = create_tree(A22, (n+1, n+1), r, ϵ) # ask @integraledelebesgue why @view does not work here
    return node
end

function count_total_tree_error(node)
    if node.is_zeros || min(size(node.matrix)...) < 16
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
    # println(size ^ 2 * zeros_percent ÷ 100)
    for i in 1:(size ^ 2 * zeros_percent ÷ 100)
        index = rand(1:size, 2)
        while matrix[index[1], index[2]] == 0
            index = rand(1:size, 2)
        end
        matrix[index[1], index[2]] = 0
    end
    return matrix
end

function compare_matrixes(mat1, mat2)
    return sum((mat1 - mat2) .^ 2)    
end

end # SVDTreeCompression