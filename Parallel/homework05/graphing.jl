using Plots



CPU = [.00467, .00467, .00467, .00467, .00467]

log2x = [ 4, 5, 6, 7, 8 ]
ELL = [.004565, 0.004325,  0.002852, 0.00242, 0.004000]
CSR = [0.002924, 0.003159, 0.002383, 0.001765, 0.003550]

plot(title = "Averge runtimes of spmv kernels by number of threads", ylabel = "time in seconds", xlabel = "log base 2 # of threads")

plot!(log2x, CPU, label = "CPU")
plot!(log2x, ELL, label = "ELL")
plot!(log2x, CSR, label = "CSR")

png("thread_runtime")

