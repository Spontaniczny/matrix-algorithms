module Determinant

using Decomposition: lu_decompose

function determinant(matrix::Matrix)::Float64
    L, U = lu_decompose(matrix)
    det = 1
    for i in 1:size(matrix, 1)
        det = det * L[i, i] * U[i, i]
    end
    return det
end

end # Determinant