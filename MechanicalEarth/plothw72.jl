using Plots

a = 1:1:1000
P = 1
r1 = .023
r3 = .026


k_1 = sqrt.(a) .* (P .+ a./pi .* (r1 - r3)) 
                             
plt1 = plot(title = "emplaced dike", xlabel = "a (meters)", ylabel="stress intensity factor", legend=false)
plot!(a, k_1)

png("hw7plot3")
