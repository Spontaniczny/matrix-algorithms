module Speed
export benchmark, Result

using LinearAlgebra: det
using TSVD
using Base.Threads
using GFlops
using SVDTreeCompression: create_tree, get_random_nonzero_matrix

struct Result
    header::Vector{Symbol}
    data::Matrix{Float64}
end

function hstack(vectors::Channel{Vector{Float64}})::Matrix{Float64}
    reduce(hcat, vectors)
end

function cases(sizes::AbstractArray{Int, 1}, percentages::AbstractArray{Int, 1})::Channel{Tuple{Int, Int, Matrix{Float64}}}
    Channel{Tuple{Int, Int, Matrix{Float64}}}() do channel
        for size in sizes
            for percentage in percentages
                put!(channel, (size, percentage, get_random_nonzero_matrix(size, percentage)))
            end
        end
    end
end

const headers::Dict{Symbol, Vector{Symbol}} = Dict(
    :time => [:size, :function, :time],
    :flops => [:size, :function, :add, :mul]
)

const additions = [
    :add16, :add32, :add64, :sub16, :sub32, :sub64
]

const multiplications = [
    :mul16, :mul32, :mul64, :div16, :div32, :div64
]

function adds(counter::GFlops.Counter)::Int
    getfield.([counter], additions) |> sum
end

function muls(counter::GFlops.Counter)::Int
    getfield.([counter], multiplications) |> sum
end


function bench_case(matrix::Matrix, percentage::Int64, rank::Int64, delta::Float64)
    time = @elapsed create_tree(matrix, (1, 1), rank, delta)
    return [size(matrix, 1), percentage, rank, delta, time]
end

function benchmark(sizes::AbstractArray{Int, 1}, zero_percentages::AbstractArray{Int, 1}, n_evals::Int)::Result
    data = Channel{Vector{Float64}}() do results
        for (n, p, matrix) in cases(sizes, zero_percentages)
            for i in 1:n_evals
                println(n, ", ",  p)
                # time = @elapsed create_tree(matrix)
                # put!(results, [size(data, 1), i, time])
                U, d, V = tsvd(matrix, size(matrix)[begin] - 1)
                put!(results, bench_case(matrix, p, 1, d[1]))
            end
        end
    end
    Result(
        data |> hstack |> transpose
    )
end

end# module