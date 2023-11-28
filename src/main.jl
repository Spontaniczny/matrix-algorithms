push!(LOAD_PATH, @__DIR__)
using Imports
Imports.@load_src_directories(".")

using Determinant: determinant
using Decomposition: lu_decompose
using Inverse: inverse

using Speed: benchmark, Result
using Base.Threads
using DataFrames: DataFrame
using CSV

function to_dataframe(result::Result)::DataFrame
    DataFrame(result.data, result.header)
end 

const path = "data/flops.csv"

function save(df::DataFrame)
    open(path, read=true, truncate=true) do io
        CSV.write(io, df)
    end
end

function main()
    functions = [lu_decompose, determinant, inverse]
    domain = 2 .^ collect(2:2)

    benchmark(functions, domain, 1, :flops) |> 
        to_dataframe |> 
        save
end

main()
