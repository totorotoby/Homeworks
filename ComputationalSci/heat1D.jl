#heat1D.jl
#=

Solves:

uₜ = kuₓₓ + F(x,t)    0 ≤ x ≤ L, 0 ≤ t ≤ T
u(x,0) = f(x)
u(0,t) = u(L,t) = 0

Using standard centered difference in space (second and fourth order) and both forward and backwards euler in time.

=#


using LinearAlgebra
using BandedMatrices
using Plots
using Printf
include("timesteppers.jl")

let

    ### parameters ###

    # solving in [0,L]
    L = 1.0
    # solving in [0,T]
    T = 2.0
    # stability parameter for timestep and spatial step
    #λ = .1

    iters = 1
    error₁ = Array{Float64,1}(undef, iters)
    error₂ = Array{Float64,1}(undef, iters)
    time₁ = Array{Float64,1}(undef, iters)
    time₂ = Array{Float64,1}(undef, iters)
    mem₁ = Array{Float64,1}(undef, iters)
    mem₂ = Array{Float64,1}(undef, iters)
    
    Δx = .1
    # iterations
    for i = 1:iters
        
        #spatial step
        #Δx = Δx/2
        @show Δx
        # spatial domain
        X = 0:Δx:L
        # spatial nodes
        n = length(X)
        @show n
        # time step
        Δt = Δx^2 * λ
        @show Δt
        # time domain
        t = 0:Δt:T
        #number of timesteps
        M = length(t)
        @show M
        # thermal diffusivity
        k = 2

        ### spatial discretization ###

        # analytic solution
        uₑ(x,t) = ℯ^(-k*t) * sin(π*x)L
        # intial condition
        f(x) = sin(π*x)
        # source function
        F(x,t) = 2 * (π^2 - 1) * ℯ^(-k*t) * sin(π*x)
        # source vector
        c(t) = [F(x,t) for x in X[2:n-1]]
        
        # stencil
        # 2nd order
        A₂ = BandedMatrix{Float64}(undef, (n-2,n-2), (1,1))
        A₂[band(0)] .= -2/Δx^2
        A₂[band(1)] .= A₂[band(-1)] .= 1/Δx^2
        A₂ = k .* A₂
        
        # solution vector (starting with intial condition)
        u = [f(x) for x in X]
        
        ### time discretization ###
        
        # make sure compiled first
        eulerF(A₂, c, u, Δt, M)
        eulerB(A₂, c, u, Δt, M)
        
        # integrate system explicitly
        R₁ = @timed eulerF(A₂, c, u, Δt, M)
        U₁ = R₁[1]
        time₁[i] = R₁[2]
        mem₁[i] = R₁[3]
        # integrate system implicitly
        R₂ = @timed eulerB(A₂, c, u, Δt, M)
        U₂ = R₂[1]
        time₂[i] = R₂[2]
        mem₂[i] = R₂[3]
        
        error₁[i] = √Δx * norm(U₁[:,end] - uₑ.(X,T))
        error₂[i] = √Δx * norm(U₂[:,end] - uₑ.(X,T))

        if i > 1
            ratio₁ = error₁[i-1]/error₁[i]
            rate₁ = log(2, ratio₁)
            ratio₂ = error₂[i-1]/error₂[i]
            rate₂ = log(2, ratio₂)
            
            @printf("convergence rate foward euler: %f\n",  rate₁)
            @printf("convergence rate backward euler: %f\n",  rate₂)
            
        end


        display(U₁)
        
        ### ploting solutions ###
        #max_y = maximum(U₁)
        #max_plot = max_y + max_y/10
        xfine = 0:.001:L

        #for m = 1:M
        plot(title="final time with Δt=.001")
        plot!(X,U₁[:,end], #=ylim = (0, max_plot)=#)
        plot!(xfine, uₑ.(xfine, T#=(m-1)*Δt=#))
        sleep(1)
        png("time3") 
        #end

    end


    ### plotting timing and mem ###
    #=
    plot(title="Times for timestepping loop with varying # spatial nodes", xlabel = "-Log₂ Δx", ylabel = "Log₂ time (seconds)")
    
    plot!(collect(1:iters), log.(2,time₁), label = "Forward Euler")
    plot!(collect(1:iters), log.(2,time₂), label = "Backwards Euler")  
    png("times_plot")

    plot(title="Memory allocated for timestepping loop with varying # spatial nodes", xlabel = "-Log₂ Δx", ylabel ="Bytes allocated")
    
    plot!(collect(1:iters), mem₁, label = "Forward Euler")
    plot!(collect(1:iters), mem₂, label = "Backwards Euler")  

    png("mem_plot")
    =#


end
