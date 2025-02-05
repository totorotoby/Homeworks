using Plots
using DelimitedFiles

times = readdlm("test_threads", ',')'
display(times)
#times = log.(2, times)

num_samples = 2:28
labels = ["serial","O(Nlog(N)) implementation","O(N) implementation"]

plot(legend=:topleft, title = "Avg over 10 Runtimes by # threads on 28 CPUs", xlabel="# of threads", ylabel=" Runtime (seconds)")

for i = 1:size(times)[1]

    plot!(num_samples, times[i, :], label=labels[i])
   
end

#gui()
png("hw2_threads")

