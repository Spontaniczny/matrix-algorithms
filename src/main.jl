push!(LOAD_PATH, @__DIR__)
using Imports
Imports.@load_src_directories(".")

using Determinant: determinant
using LinearAlgebra: lu


function main()
    println(determinant(A))
end

main()