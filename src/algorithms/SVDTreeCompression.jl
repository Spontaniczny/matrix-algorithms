module SVDTreeCompression

using TSVD
using LinearAlgebra: Diagonal, svd
using Colors: Gray
using Plots: plot, savefig

export create_tree, count_absolut_tree_error, count_relative_tree_error, draw_tree, compare_matrixes

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
        return node
    end

    # julia's tsvd is broken, so we use normal svd, but we take min s values < ϵ
    # A += Diagonal(repeat([1e-15],  min(size(node.matrix)...)))
    # u, s, v = tsvd(A, r + 1)
    # A -= Diagonal(repeat([1e-15],  min(size(node.matrix)...)))

    u, s, v = svd(A)

    # if s[r + 1] < ϵ || min(size(node.matrix)...) <= 4 to jest jakby force compress. Chcemy tak robić?
    i = min(r + 1, size(s)[1])
    if s[i] < ϵ
        # to use the smallest value < ϵ
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

function count_absolut_tree_error(node)
    if node.is_zeros || min(size(node.matrix)...) <= 2
        return 0
    elseif node.is_compressed
        return compare_matrixes(node.matrix, node.u * Diagonal(node.s) * transpose(node.v))
    end

    error = 0
    error += count_absolut_tree_error(node.top_right_child)
    error += count_absolut_tree_error(node.top_left_child)
    error += count_absolut_tree_error(node.bottom_left_child)
    error += count_absolut_tree_error(node.bottom_right_child)
    return error
end

function count_relative_tree_error(node, abs_error)
    return abs_error / (size(node.matrix)[1] ^ 2)
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
            sp = node.top_left #start point
            node_size = size(node.matrix)[1]
            matrix[sp[1]:sp[1]+node_size-1, sp[2]:sp[2]+node_size-1] .= 0
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
    # plot(Gray.(matrix)) # it works only in console afaik
    return matrix
end

# matrix = get_random_nonzero_matrix(128, 99)
# s = svd(matrix)
# for i in 1:128
#     u, s, v = tsvd(matrix, i)
#     # u, s, v = svd(matrix)
#     # matrix2 = u[:, 1:i] * Diagonal(asd[1:i]) * transpose(v[:, 1:i])
#     matrix2 = u * Diagonal(s) * transpose(v)
#     # println(i, " ", s[i], " ", asd[i])
#     println(i, " ", compare_matrixes(matrix, matrix2))
#     println()
# end


# m_size_tab = [256, 512, 1024]
# zeros_percent_tab = [99, 98, 95, 90, 80] 

# _, s, _ = svd(matrix)

# s_tab = [s[2], s[length(s) ÷ 2], s[end]]
# b_tab = [1, 4]

# s = s[2]
# s = s[length(s) ÷ 2]
# s = s[end]

# for m_size in m_size_tab
#     for zeros_percent in zeros_percent_tab
#         for b in b_tab
#             for s in s_tab
#                 println(m_size, zeros_percent, b, s)
#                 matrix = get_random_nonzero_matrix(m_size, zeros_percent)
#                 xd = create_tree(matrix, (1, 1), b, s)
#                 # x = count_absolut_tree_error(xd)
#                 # x_rel = count_relative_tree_error(xd, x)
#                 xd2 = draw_tree(xd)
#                 # println(x)
#                 # println(x_rel)
#                 plot(Gray.(xd2), title = "$(m_size)x$(m_size), $(zeros_percent)%, b=$b, ϵ=$(round(s, digits = 5))", dpi = 600, grid = false)
#                 savefig("benchmarks/lab3/$(m_size)x$(m_size)_$(zeros_percent)%_b=$(b)_eps=$(round(s, digits = 5)).png")
#             end
#         end
#     end
# end

end # SVDTreeCompression