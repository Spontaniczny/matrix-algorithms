push!(LOAD_PATH, @__DIR__)
using Imports
Imports.@load_src_directories(".")

function main()
    print("Hello World!")
end

main() |> display