module Helper
export generate_spd_matrix, is_correct, measure_time, results_in_table, make_plot, analyze_results

using LinearAlgebra, Statistics, Plots, Printf, Base.Threads, Measures

function generate_spd_matrix(n)
    A = rand(n, n)
    A = A * A' + I(n)
    return A
end

function is_correct(A, L)
    return norm(A - L * L') < 1e-8
end

function measure_time(f, A; iterations=3)
    times = Float64[]
    for _ in 1:iterations
        t_start = time()
        f(A)
        t_end = time()
        push!(times, t_end - t_start)
    end
    
    return minimum(times)
end

function results_in_table(results, speedup, sizes)
    println("РЕЗУЛЬТАТЫ ЭКСПЕРИМЕНТОВ")
    println("="^80)
    println("\nСреднее время выполнения (секунды):")
    println("\n" * "-"^80)
    println(@sprintf("%-10s | %-12s | %-12s | %-12s | %-12s | %-12s", 
                     "Размер", "Serial", "ThOuter", "ThInner", "PolOuter", "PolInner"))
    println("-"^80)
    
    for n in sizes
        t_serial = mean(results["serial"][n])
        t_th_outer = mean(results["threads_outer"][n])
        t_th_inner = mean(results["threads_inner"][n])
        t_pol_outer = mean(results["polyester_outer"][n])
        t_pol_inner = mean(results["polyester_inner"][n])
        
        println(@sprintf("%-10d | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f",
                         n, t_serial, t_th_outer, t_th_inner, t_pol_outer, t_pol_inner))
    end
    
    println("\n" * "-"^80)
    println("\nУскорение (по сравнению с последовательной версией):")
    println("\n" * "-"^80)
    println(@sprintf("%-10s | %-12s | %-12s | %-12s | %-12s", 
                     "Размер", "ThOuter", "ThInner", "PolOuter", "PolInner"))
    println("-"^80)
    
    for n in sizes
        s_th_outer = speedup["threads_outer"][n]
        s_th_inner = speedup["threads_inner"][n]
        s_pol_outer = speedup["polyester_outer"][n]
        s_pol_inner = speedup["polyester_inner"][n]
        
        println(@sprintf("%-10d | %-12.2f | %-12.2f | %-12.2f | %-12.2f",
                         n, s_th_outer, s_th_inner, s_pol_outer, s_pol_inner))
    end
    
    println("-"^80)
end

function analyze_results(speedup, sizes)
    println()
    println("="^80)
    println("НАИБОЛЕЕ ЭФФЕКТИВНАЯ СТРАТЕГИЯ ПАРАЛЛЕЛИЗАЦИИ")
    println("="^80)
    
    for n in sizes
        max_speedup = 0.0
        best_method = ""
        
        for method in keys(speedup)
            if haskey(speedup[method], n)
                s = speedup[method][n]
                if s > max_speedup
                    max_speedup = s
                    best_method = method
                end
            end
        end
        
        println("\nДля матрицы размера $n:\n")
        println("  Лучшая стратегия: $best_method с ускорением $(@sprintf("%.2f", max_speedup))\n")
        println("  Эффективность: $(@sprintf("%.1f", (max_speedup / nthreads()) * 100))%\n")
        println("="^80)
    end
end

# TODO: Выделить построение графиков в отдельные методы
function make_plots(results, speedup, sizes)
    println("ПОСТРОЕНИЕ ГРАФИКОВ")
    println("="^80)
    p1 = plot(
        title="Время выполнения метода Холецкого",
        xlabel="Размер матрицы (n)",
        ylabel="Время (секунды)",
        xscale=:identity,
        yscale=:log10,
        legend=:bottomright,
        grid=true,
        minorgrid=true,
        xticks=sizes,
    )
    
    avg_times = Dict()
    for (method, data) in results
        avg_times[method] = [mean(data[n]) for n in sizes]
    end
    
    colors = palette(:tab10)

    plot!(p1, sizes, avg_times["serial"], 
          label="Последовательная", 
          marker=:circle, linewidth=2, color=colors[1])
    plot!(p1, sizes, avg_times["threads_outer"], 
          label="Threads (внешний)", 
          marker=:square, linewidth=2, color=colors[2])
    plot!(p1, sizes, avg_times["threads_inner"], 
          label="Threads (внутренний)", 
          marker=:diamond, linewidth=2, color=colors[3])
    plot!(p1, sizes, avg_times["polyester_outer"], 
          label="Polyester (внешний)", 
          marker=:utriangle, linewidth=2, color=colors[4])
    plot!(p1, sizes, avg_times["polyester_inner"], 
          label="Polyester (внутренний)", 
          marker=:dtriangle, linewidth=2, color=colors[5])

    # min_speedup = minimum([minimum([speedup[m][n] for n in sizes]) for m in keys(speedup)])
    # max_speedup = maximum([maximum([speedup[m][n] for n in sizes]) for m in keys(speedup)])

    p2 = plot(
        title="Ускорение параллельных версий",
        xlabel="Размер матрицы (n)",
        ylabel="Ускорение (раз)",
        xscale=:identity,
        xticks=sizes,
        legend=:topright,
        grid=true,
        ylim=(0, nthreads()),
        xlim=(minimum(sizes)-2, maximum(sizes)+2)
    )
    
    hline!([nthreads()], label="Идеальное ($(nthreads()) потока)", 
           linestyle=:dash, color=:gray, linewidth=1)

    plot!(p2, sizes, [speedup["threads_outer"][n] for n in sizes], 
          label="Threads (внешний)", marker=:square, linewidth=2, color=colors[2])
    plot!(p2, sizes, [speedup["threads_inner"][n] for n in sizes], 
          label="Threads (внутренний)", marker=:diamond, linewidth=2, color=colors[3])
    plot!(p2, sizes, [speedup["polyester_outer"][n] for n in sizes], 
          label="Polyester (внешний)", marker=:utriangle, linewidth=2, color=colors[4])
    plot!(p2, sizes, [speedup["polyester_inner"][n] for n in sizes], 
          label="Polyester (внутренний)", marker=:dtriangle, linewidth=2, color=colors[5])
    
    savefig(p1, "cholesky_times.png")
    println("\nГрафик сохранен в cholesky_times.png\n")
    savefig(p2, "cholesky_speedup.png")
    println("График сохранен в cholesky_speedup.png\n")

    combined_plot = plot(p1, p2, layout=(2, 1), size=(1000, 1200))
    savefig(combined_plot, "cholesky_results.png")
    println("График сохранен в cholesky_results.png\n")
    println("="^80)
    println("ЭКСПЕРИМЕНТЫ ЗАВЕРШЕНЫ")
    println("="^80)
end

end
