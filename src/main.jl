push!(LOAD_PATH, @__DIR__)
using Imports
Imports.@load_src_directories(".")

using Determinant: determinant
using Decomposition: lu_decompose
using Inverse: inverse
using SVDTreeCompression: create_tree, get_random_nonzero_matrix


using Speed: benchmark, Result
using Base.Threads
using DataFrames: DataFrame
using CSV
using LinearAlgebra: inv

function to_dataframe(result::Result)::DataFrame
    DataFrame(result.data, result.header)
end 

const path = "data/flops.csv"

function save(df::DataFrame)
    open(path, read=true, truncate=true) do io
        CSV.write(io, df)
    end
end

function force_precompilation()::Nothing
    @sync begin
        @spawn create_tree(get_random_nonzero_matrix(16, 90))
    end

    nothing
end

function main()
    force_precompilation()

    sizes = 2 .^ collect(4:12)
    zero_percentages = [80, 90, 95, 98, 99]

    benchmark(sizes, zero_percentages, 1) |> 
        to_dataframe |> 
        save
end

# Commend for @Kuba because he always forget what to type in console
# julia --project=. src/main.jl

main()
