module Plot
    using Revise
    using MatrixAlgorithms
    using CSV, Tables

    function compute_total_operations(sizes, algorithm, samples=3)
        operations = []

        for n in sizes
            additions = 0
            multiplications = 0
            for sample in 1:samples
                A, B = AugmentMatrix(Float64.(rand(n, n))), AugmentMatrix(Float64.(rand(n, n)))
                solution = algorithm(A, B)
                additions += solution.add
                multiplications += solution.mul
            end
            
            data_point = (size=n, add=additions/samples, mul=multiplications/samples)
            push!(operations, data_point)
        end
        return operations
    end

    function operations_to_csv(operations, filepath)
        sizes = map(op -> op.size, operations)
        adds = map(op -> op.add, operations)
        muls = map(op -> op.mul, operations)
        
        for data_point in operations
            CSV.write(filepath, (size = sizes, additions=adds, multiplications=muls))
        end
    end
end # module Plot