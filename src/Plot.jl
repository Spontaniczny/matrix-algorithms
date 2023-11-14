module Plot
    using Revise
    using CSV, Tables
    using BenchmarkTools
    using MatrixAlgorithms

    export compute_total_operations, operations_to_csv, compute_total_times, times_to_csv

    function compute_total_operations(sizes, algorithm, samples=3, verbose=false)
        operations = []

        for n in sizes
            additions = 0
            multiplications = 0

            if verbose
                print("Wokring on n = ", n)
            end

            for sample in 1:samples
                A, B = AugmentMatrix(Float64.(rand(n, n))), AugmentMatrix(Float64.(rand(n, n)))
                solution = algorithm(A, B)
                additions += solution.add
                multiplications += solution.mul

                if verbose 
                    print("...", sample)
                end
            end
            
            if(verbose)
                println("")
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


    function compute_total_times(sizes, algorithm, verbose=false)
        times = []

        for n in sizes

            if verbose
                print("Wokring on n = ", n)
            end

            
            A, B = AugmentMatrix(Float64.(rand(n, n))), AugmentMatrix(Float64.(rand(n, n)))
            compute_time = round(@benchmark algorithm(A, B) / 10e5)
            
            
            if(verbose)
                println("")
            end            

            data_point = (size=n, cumpute_time=compute_time)
            push!(times, data_point)
        end
        return times
    end

    function times_to_csv(times, filepath)
        sizes = map(op -> op.size, times)
        op_times = map(op -> op.op_time, times)
        
        for data_point in operations
            CSV.write(filepath, (size = sizes, time = op_times))
        end
    end
end # module Plot