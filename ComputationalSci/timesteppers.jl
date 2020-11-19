using LinearAlgebra


function eulerF(A, c::Function, u₀::Array{Float64,1},
                Δt::Float64, M::Integer)

    N = length(u₀)
    # solution matrix
    U = zeros(N, M)
    U[:, 1] .= u₀[:]
    
    for m = 2:M
        # take an euler step
        # (previous solution + timestep * change in solution)
        U[2 : N-1, m] = U[2 : N-1, m-1] .+ Δt*(A*U[2: N-1, m-1] .+ c((m)*Δt))
        @show m
        @show U[3, m]
    end
    
    return U

end


function eulerB(A, c::Function, u₀::Array{Float64,1},
                Δt::Float64, M::Integer)

    N = length(u₀)
    # solution matrix
    U = zeros(N,M)
    U[:,1] .= u₀[:]

    for m = 2:M
        B = (I - Δt .* A)
        F = lu(B)
        b = U[2 : N-1, m-1] + Δt * c((m)*Δt)
        U[2 : N-1, m] = F\b
    end
    
    return U
end
    
