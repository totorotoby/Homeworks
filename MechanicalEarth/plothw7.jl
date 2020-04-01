using Plots

theta =-pi : .01 : pi

sigma_rr = -sin.(theta./2) .* (1 .- 3 .*sin.(theta./2) .^ 2)
sigma_tt = 3 .* sin.(theta./2) .* cos.(theta./2) .^ 2
sigma_rt = -cos.(theta./2) .* (1 .- 3sin.(theta./2) .^ 2)
                             
plt1 = plot(title = "Mode 2 Stress Components", xlabel = "theta", ylabel="stress")
plot!(theta, sigma_rr, label="sigma rr")
plot!(theta, sigma_tt, label="sigma tt")
plot!(theta, sigma_rt, label="sigma rt")

png("hw7plot2")
