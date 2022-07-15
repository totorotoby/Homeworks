using Plots
using DelimitedFiles

times = readdlm("num_threads", ',')

@show length(times)


threads = collect(1:length(times))


plot(title="runtime of quadtrature by number of threads", xlabel = "#number of threads", ylabel = "runtime (seconds)")
plot!(threads, times)

png("threads_runtime")
