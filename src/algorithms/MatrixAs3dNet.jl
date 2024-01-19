module MatrixAs3dNet

using LinearAlgebra: Diagonal
using Colors: Gray
using Plots: plot, savefig
using DataStructures: Deque

using ..SVDTreeCompression: create_tree, count_absolut_tree_error, count_relative_tree_error, draw_tree, compare_matrixes




function get_random_3d_net_matrix(k)
    n = 2^k
    n_2 = n^2 

    matrix = zeros(n^3, n^3)
    matrix_size = size(matrix)[1]

    for vertex in 1:matrix_size
        level = (vertex - 1) ÷ n_2
        rest = (vertex - 1) % n_2
        row = rest ÷ n
        col = rest % n

        matrix[vertex, vertex] = rand()
        if level > 0
            top_level_neighbor = vertex - n_2
            matrix[vertex, top_level_neighbor] = rand()
        end

        if level < n - 1
            bottom_level_neighbor = vertex + n_2
            matrix[vertex, bottom_level_neighbor] = rand()
        end

        if row > 0
            top_neighbor = vertex - n
            matrix[vertex, top_neighbor] = rand()
        end

        if row < n - 1
            bottom_neighbor = vertex + n
            matrix[vertex, bottom_neighbor] = rand()
        end

        if col > 0
            left_neighbor = vertex - 1
            matrix[vertex, left_neighbor] = rand()
        end

        if col < n - 1
            right_neighbor = vertex + 1
            matrix[vertex, right_neighbor] = rand()
        end
    end

    return matrix
end

function get_adjacency(matrix)
    s = size(matrix)[1]

    adjency = Dict{Int, Set}()
    for i in 1:s
        adjency[i] = Set()
    end

    for i in 1:s
        for j in 1:s
            if i != j && matrix[i, j] > 0
                push!(adjency[i], j)
            end
        end
    end
    return adjency
end

function min_degree(matrix)
    s = size(matrix)[1]

    result = []
    adjacency = get_adjacency(matrix)

    for i in 1:s
        best_vertex = -1
        best_value = s + 1

        for (vertex, neighbors) in adjacency
            if length(neighbors) < best_value
                best_value = length(neighbors)
                best_vertex = vertex
            end
        end

        for vertex in keys(adjacency)
            delete!(adjacency[vertex], best_vertex)
        end

        for neighbor in adjacency[best_vertex]
            union!(adjacency[neighbor], setdiff(adjacency[best_vertex], Set(neighbor)))
        end

        delete!(adjacency, best_vertex)
        push!(result, best_vertex)
    end

    return result
end

function cuthill_mckee(matrix)

    function bfs()
        while length(q) != 0
            vertex = popfirst!(q)

            if visited[vertex]
                continue
            end

            visited[vertex] = true
            push!(result, vertex)

            neighbors = sort(collect(adjacency[vertex]), by=length)

            for neighbor in neighbors
                if !visited[neighbor]
                    push!(q, neighbor)
                end
            end
        end
    end
        
    s = size(matrix)[1]

    adjacency = get_adjacency(matrix)
    degrees = []
    for (vertex, neighbors) in adjacency
        push!(degrees, (vertex, length(neighbors)))
    end

    sorted_vertices = sort(degrees, by = x -> x[2])

    for i in 1:length(sorted_vertices)
        sorted_vertices[i] = sorted_vertices[i][1]
    end

    result = []

    q = Deque{Int}()
    visited = repeat([false], s)

    for vertex in sorted_vertices
        if !visited[vertex]
            push!(q, vertex)
            bfs()
        end
    end
    return result
end

function reversed_cuthill_mckee(matrix)
    return reverse(cuthill_mckee(matrix))  # maybe use reverse!
end


function apply_permutation(matrix, permutation)
    result = copy(matrix)
    
    for i in 1:length(permutation)
        if i != permutation[i]
            result[i, :] = matrix[permutation[i], :]
        end
    end

    matrix = copy(result)

    for i in 1:length(permutation)
        if i != permutation[i]
            result[:, i] = matrix[:, permutation[i]]
        end
    end

    return result
end

# m = get_random_3d_net_matrix(4)
# perm_min_degree = reversed_cuthill_mckee(m)
# matrix_min_degree = apply_permutation(m, perm_min_degree)
# xd = create_tree(matrix_min_degree, (1, 1), 1, 0.1)
# xd2 = draw_tree(xd)

# varinfo
end