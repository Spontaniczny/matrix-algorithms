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
        is_zeros = false
        is_compressed = false
        return foo
    end
end

function create_tree(A::MatrixOrView, top_left::Tuple{Int, Int} = (1, 1), r::Int64 = 1, ϵ::Float64 = 1.0)
    node = Tree_Node(A, top_left)
    node.is_zeros = iszero(A)
    if node.is_zeros || min(size(node.matrix)...) <= 2
    # if node.is_zeros
    #     return node
    # elseif min(size(node.matrix)...) < 16
        return node
    end

    # is_ok = false
    # counter = 0
    
    # while !is_ok && counter < 100
    #     local u, s, v # or define to nothing
    #     try
    #         u, s, v = tsvd(A, r + 1)
    #         is_ok = true
    #     catch e
    #         x, _ = size(A)
    #         index = rand(1:x, 2)
    #         A[index[1], index[2]] += 1e-15
    #         counter += 1
    #     end
    # end

    A += Diagonal(repeat([1e-15],  min(size(node.matrix)...)))
    u, s, v = tsvd(A, r + 1)
    A -= Diagonal(repeat([1e-15],  min(size(node.matrix)...)))
    # if s[r + 1] < ϵ || min(size(node.matrix)...) <= 4 to jest jakby force compress. Chcemy tak robić?
    if s[r + 1] < ϵ
        i = r + 1
        while s[i] < ϵ && i != 1
            i -= 1
        end
        node.is_compressed = true
        node.u = u[:, 1:i]
        node.s = s[1:i]
        node.v = v[:, 1:i]  # na pseudokodzie jest tu mnozenie Diagonal(node.s) * transpose(node.v)
        return node
    end
    A11, A12, A21, A22 = split_view(A)
    n, _ = div.(size(A), 2)

    node.is_compressed = false

    node.top_left_child = create_tree(A11, node.top_left, r, ϵ) # ask @integraledelebesgue why @view does not work here
    node.top_right_child = create_tree(A12, (node.top_left[1], node.top_left[2] + n), r, ϵ) # ask @integraledelebesgue why @view does not work here
    node.bottom_left_child = create_tree(A21, (node.top_left[1] + n, node.top_left[2]), r, ϵ) # ask @integraledelebesgue why @view does not work here
    node.bottom_right_child = create_tree(A22, (node.top_left[1] + n, node.top_left[2] + n), r, ϵ) # ask @integraledelebesgue why @view does not work here
    return node
end

function count_total_tree_error(node)
    if node.is_zeros || min(size(node.matrix)...) <= 2
        return 0
    elseif node.is_compressed
        return compare_matrixes(node.matrix, node.u * Diagonal(node.s) * transpose(node.v))
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

function draw_tree(tree)

    function fill_matrix(matrix, node)

        if node.is_zeros
            # println("!@#")
            return
        elseif min(size(node.matrix)...) <= 2
            # println("XD")
            return
        elseif node.is_compressed
            sp = node.top_left #start point
            node_size = size(node.matrix)[1]
            no_ranks = size(node.s)[1]
            # println(sp)
            # println(node_size)

            for i in 0:no_ranks-1
                matrix[sp[1]:sp[1]+node_size-1, sp[2]+i] .= 0
                matrix[sp[1]+i, sp[2]:sp[2]+node_size-1] .= 0
            end

            return
        end
        
        fill_matrix(matrix, node.top_right_child)
        fill_matrix(matrix, node.top_left_child)
        fill_matrix(matrix, node.bottom_left_child)
        fill_matrix(matrix, node.bottom_right_child)
        return
    end

    matrix = ones(size(tree.matrix))
    fill_matrix(matrix, tree)
    plot(Gray.(matrix))
    return matrix
end

matrix = get_random_nonzero_matrix(256, 98)
# for i in 1:128
#     u, s, v = tsvd(matrix, i)
#     matrix2 = u * Diagonal(s) * transpose(v)
#     println(i, " ", s[i])
#     println(i, " ", compare_matrixes(matrix, matrix2))
#     println()
# end



# xd = create_tree(matrix, (1, 1), 4, 0.001)
# x = count_total_tree_error(xd)
# xd2 = draw_tree(xd)

# println(x)

end # SVDTreeCompression