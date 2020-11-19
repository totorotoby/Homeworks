using SparseArrays
using LinearAlgebra
using IterativeSolvers
using Plots
include("my_solvers.jl")

let
    
    N = 1000

    # dense implementation
    A₁= zeros(N,N)
    for i = 1:N
        A₁[i,i] = -2.0
    end
    for i = 2:N
        A₁[i,i-1] = 1.0
        A₁[i-1, i] = 1.0
    end
    
    A₁copy = deepcopy(A₁)

    # sparse implemenation
    i1 = [i for i=1:N]
    i2 = [i for i=2:N]
    i3 = [i for i=1:N-1]
    v1 = [-2.0 for v=1:N]
    v2 = [1.0 for v=1:N-1]
    I = vcat(i1, i2, i3)
    J = vcat(i1, i3, i2)
    V = vcat(v1, v2, v2)
    
    A₂ = sparse(I,J,V)
    A₂copy = deepcopy(A₂)

    x₀ = zeros(N)
    
    b = rand(N)
    
    lu(A₁)

    println("Julia LU factorization on dense matrix:")
    @time F = lu(A₁)
    println("\njulia LU solver on dense matrix:")
    @time x = F\b

    
    println("\nMy LU factorization on dense matrix:")
    @time P = computeLUP!(A₁copy)
    println("\nMy LU solver on dense matrix:")
    @time x = LUPsolve(A₁copy, P, b)
    

    println("\nMy conjugate gradient on dense matrix:")
    @time x = conj_grad(A₁, b, x₀, 1e-15, 2000)
    
    
    println("\nJulia LU factorization on sparse matrix:")
    @time F = lu(A₂)
    println("\njulia LU solver  on sparse matrix:")
    @time x = F\b

    
    println("\nMy LU factorization on sparse matrix:")
    @time P = computeLUP!(A₂copy)
    println("\nMy LU solver  on sparse matrix:")
    @time x = LUPsolve(A₁copy, P, b)
    

    println("\nMy conjugate gradient on sparse matrix:")
    @time x = conj_grad(A₂, b, x₀, 1e-15, 2000)

    println("\n\n_____________________________________\n\n")
    
    
    for N = 10 * 10 .^ (0:2)
    
        A = rand(N,N)
        # symmetrify
        A = (1/2)*(A+A')
        # positivedefify
        A = A'*A
        cx = zeros(N)
        b = rand(N)

        
        println("N:", N)
        println("LU factorization:")
        @time F = lu(A)
        println("LU solve:")
        @time x = F\b
        println("CG solve:")
        @time cg!(cx, A, b; tol=1e-8)
        println("My CG solve:")
        @time conj_grad(A, b, cx, 1e-8, 2*N)
        println("\n\n\n")
        
    end
    
    N = 1000

    A = rand(N,N)
    # symmetrify
    A = (1/2)*(A+A')
    # positivedefify
    A = A'*A
    cx = zeros(N)
    b = rand(N)

    ϵ = 1 * 10.0 .^ (-1:-1:-12)
    ϵ = collect(ϵ)
    times = zeros(length(ϵ))
    x_axis = [i for i=1:12]
    @show ϵ
    for (i,e) = enumerate(ϵ)
        time = @elapsed conj_grad(A, b, cx, e, 10*N)
        times[i] = time
    end
  
    plot(x_axis, times, xlabel="Log₁₀ error tolerance", ylabel="time (seconds)", title="Time to Run conj_grad with Decreasing Relative Error")
    png("error_plot")
    
    
    nothing
    




end
