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
    a = rand(8, 8)
    @sync begin
        @spawn lu_decompose(a)
        @spawn determinant(a)
        @spawn inverse(a)
    end

    nothing
end

function main()
    # force_precompilation()

    functions = [lu_decompose, determinant, inv]
    domain = 2 .^ collect(1:9)

    benchmark(functions, domain, 1, :flops) |> 
        to_dataframe |> 
        save
end

main()
