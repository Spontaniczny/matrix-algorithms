push!(LOAD_PATH, @__DIR__)
using Imports
Imports.@load_src_directories(".")

using Inverse: inverse!


function main()
    A = rand(2, 2)
    A_inv = inverse!(A)
    println("Real solution")
    display(inv(A))
    println("Our solution")
    display(A_inv)
end

main()