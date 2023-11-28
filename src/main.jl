push!(LOAD_PATH, @__DIR__)
using Imports
Imports.@load_src_directories(".")

using Multiplication: Strassen

function main()
    A, B = rand(2, 2), rand(2, 2)
    C = Strassen(A, B)
    println(C)
end

main()