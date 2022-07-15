using Plots
using DelimitedFiles


thread_data = readdlm("thread_data")

display(thread_data)
plot(title = "Time of CG times with 13000 iters by # threads per block", xlabel = "log2 # threads", ylabel = "time(s)", legend=false)
plot!(log.(2,thread_data[:,1]), thread_data[:,2])
png("threads")
