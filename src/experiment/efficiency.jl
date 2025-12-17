using LinearAlgebra, Statistics, Printf, Base.Threads, Plots

include(joinpath(@__DIR__, "..", "serial_version", "cholesky_serial.jl"))
include(joinpath(@__DIR__, "..", "parallel_versions", "cholesky_polyester.jl"))
include(joinpath(@__DIR__, "..", "parallel_versions", "cholesky_threads.jl"))
include(joinpath(@__DIR__, "helper.jl"))

using Main.CholeskySerial: cholesky_serial
using Main.CholeskyPolyester: cholesky_polyester_inner_safe, cholesky_polyester_outer_safe
using Main.CholeskyThreads: cholesky_threads_inner_safe, cholesky_threads_outer_safe
using Main.Helper: generate_spd_matrix, is_correct, measure_time, results_in_table, make_plots, analyze_results



function test_efficiency()
    sizes = [100]

    results = Dict(
        "serial" => Dict{Int, Vector{Float64}}(),
        "threads_outer" => Dict{Int, Vector{Float64}}(),
        "threads_inner" => Dict{Int, Vector{Float64}}(),
        "polyester_outer" => Dict{Int, Vector{Float64}}(),
        "polyester_inner" => Dict{Int, Vector{Float64}}()
    )
    
    speedup = Dict(
        "threads_outer" => Dict{Int, Float64}(),
        "threads_inner" => Dict{Int, Float64}(),
        "polyester_outer" => Dict{Int, Float64}(),
        "polyester_inner" => Dict{Int, Float64}(),
    )
    
    println("="^80)
    println("ЗАПУСК ЭКСПЕРИМЕНТОВ")
    println("="^80)
    println("ДОСТУПНО ПОТОКОВ: $(nthreads())")
    println("="^80)
    
    for (idx, n) in enumerate(sizes)
        println("\nТестирование матрицы размера $n ($idx/$(length(sizes)))\n")
        
        A = generate_spd_matrix(n)
        
        times_serial = Float64[]
        times_threads_outer = Float64[]
        times_threads_inner = Float64[]
        times_polyester_outer = Float64[]
        times_polyester_inner = Float64[]
        
        println("Проверка корректности...\n")
        
        L_serial = cholesky_serial(A)
        if !is_correct(A, L_serial)
            error("Ошибка в последовательной версии для n=$n\n")
        end
        
        versions = [
            ("Threads (внешний)", cholesky_threads_outer_safe),
            ("Threads (внутренний)", cholesky_threads_inner_safe),
            ("Polyester (внешний)", cholesky_polyester_outer_safe),
            ("Polyester (внутренний)", cholesky_polyester_inner_safe)
        ]
        
        for (name, func) in versions
            L_test = func(A)
            if !is_correct(A, L_test)
                @warn "Возможная ошибка в $name для n=$n\n"
            end
        end
        
        println("Измерение производительности...\n")
        
        for _ in 1:3
            push!(times_serial, measure_time(cholesky_serial, A))
            push!(times_threads_outer, measure_time(cholesky_threads_outer_safe, A))
            push!(times_threads_inner, measure_time(cholesky_threads_inner_safe, A))
            push!(times_polyester_outer, measure_time(cholesky_polyester_outer_safe, A))
            push!(times_polyester_inner, measure_time(cholesky_polyester_inner_safe, A))
        end
        
        results["serial"][n] = times_serial
        results["threads_outer"][n] = times_threads_outer
        results["threads_inner"][n] = times_threads_inner
        results["polyester_outer"][n] = times_polyester_outer
        results["polyester_inner"][n] = times_polyester_inner
        
        avg_serial = mean(times_serial)
        
        speedup["threads_outer"][n] = avg_serial / mean(times_threads_outer)
        speedup["threads_inner"][n] = avg_serial / mean(times_threads_inner)
        speedup["polyester_outer"][n] = avg_serial / mean(times_polyester_outer)
        speedup["polyester_inner"][n] = avg_serial / mean(times_polyester_inner)
        
        println("Среднее время выполнения:\n")
        println("  Последовательный метод: $(@sprintf("%.3f", avg_serial)) с\n")
        println("  Threads:\n")
        println("       (внешний): $(@sprintf("%.3f", mean(times_threads_outer))) с, ускорение: $(@sprintf("%.2f", speedup["threads_outer"][n]))\n")
        println("       (внутренний): $(@sprintf("%.3f", mean(times_threads_inner))) с, ускорение: $(@sprintf("%.2f", speedup["threads_inner"][n]))\n")
        println("  Polyester:\n")
        println("       (внешний): $(@sprintf("%.3f", mean(times_polyester_outer))) с, ускорение: $(@sprintf("%.2f", speedup["polyester_outer"][n]))\n")
        println("       (внутренний): $(@sprintf("%.3f", mean(times_polyester_inner))) с, ускорение: $(@sprintf("%.2f", speedup["polyester_inner"][n]))\n")
        println("="^80)
    end
    
    return results, speedup, sizes
end


results, speedup, sizes = test_efficiency()
results_in_table(results, speedup, sizes)
analyze_results(speedup, sizes)
make_plots(results, speedup, sizes)
