using TSVD
using LinearAlgebra: Diagonal
# using Decomposition: split_view


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

xd = 
matrix = get_random_nonzero_matrix(xd, 99)
n, _ = div.(size(matrix), 2)
matrix2 = @views matrix[begin:n, begin:n]
# println(matrix)

# matrix = [5 5 6 8 3; 7 4 2 3 4; 3 4 9 4 9; 4 7 9 1 5; 8 4 1 5 7]



# println(matrix)
# println(size(u))
# println(size(s))
# println(size(v))265
# println(u * Diagonal(s))
minsqrdiff = Inf
mini = 0
for i in 1:xd-1
    global minsqrdiff, mini
    u, s, v = tsvd(matrix, i)
    svdmatrix = u*Diagonal(s)*transpose(v)
    # println(round.(svdmatrix, digits = 5))
    sqrdiff = compare_matrixes(matrix, svdmatrix)
    if sqrdiff < minsqrdiff
        minsqrdiff = sqrdiff
        mini = i
    end
    println(i, " ", round.(sqrdiff, digits = 5))
end

println("min ", mini, " ", minsqrdiff)
println("3.2569094358862474 vs ", minsqrdiff * 1024)