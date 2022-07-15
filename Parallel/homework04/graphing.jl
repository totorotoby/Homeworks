using Plots
using DelimitedFiles

times = readdlm("thread_nums", ',')
@show times
display(times)


num_samples = 1:12

plot(legend=false, title = "Runtimes by # threads on 12 CPUs", xlabel="# of threads", ylabel=" Runtime (seconds)")

plot!(num_samples, times)
println("asdf")   
#gui()
png("hw4_threads")

