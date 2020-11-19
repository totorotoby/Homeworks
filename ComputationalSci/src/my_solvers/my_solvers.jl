export computeLUP, LUPsolve, conj_grad


module my_solvers

using LinearAlgebra
using Plots


"""
        computeLUP!(A)
        
        Takes in a matrix A and computes its LUP factorization storing L and U back in A,
        and returning P.
        
        ## EXAMPLE
        julia> A = rand(3,3)
        3×3 Array{Float64,2}:
         0.780188  0.819181  0.19653
         0.717729  0.643044  0.721505
         0.908919  0.667674  0.371328

        julia> computeLUP!(A)
        3×3 Array{Float64,2}:
         0.0  0.0  1.0
         1.0  0.0  0.0
         0.0  1.0  0.0

        julia> A
        3×3 Array{Float64,2}:
         0.908919  0.667674   0.371328
         0.858369  0.24607   -0.122206
         0.789651  0.470656   0.485802

        """
function computeLUP!(A)

    N = size(A)[1]
    P = Matrix(1.0I, N, N)
    #finding max value of kth column between k and nth rows
    for k=1:N
        
        # row with maximum value on kth column
        max_row = findmax(broadcast(abs, A[k:N, k]))[2] + (k-1)

        if max_row != k
            ############################################
            #  updating permutation vector by swapping #
            ############################################
            temp = P[k,:]
            P[k, :] = P[max_row,:]
            P[max_row, :] = temp
            
            ##########################################
            # row swapping matrix according to pivot #
            ##########################################
            temp = A[k, :]
            A[k, :] = A[max_row, :]
            A[max_row, :] = temp
            
        end
        
        
        # get multiplers and place in lower part of A (store L inverse instead of Lᵢ so no negative is needed)
        for i=k+1:N                
            A[i,k] /= A[k,k]
            for j=k+1:N
                A[i,j] -= A[i,k] * A[k,j]
            end
        end
    end
    return P
end


"""
        LUPsolve(A, P, b)

        Solves the system Ax=b using a LU decomposion stored in A and a permutation matrix in P.

        ## Example

        julia> b = rand(3)
        3-element Array{Float64,1}:
         0.6564446141551763
         0.039206369741312086
         0.6120722514542289

        julia> A = rand(3,3)
        3×3 Array{Float64,2}:
         0.6649     0.801754  0.300337
         0.0474491  0.719728  0.765741
         0.346125   0.922222  0.196634

        julia> P = computeLUP!(A)
        3×3 Array{Float64,2}:
         1.0  0.0  0.0
         0.0  1.0  0.0
         0.0  0.0  1.0

        julia> LUPsolve(A,P,b)
        3-element Array{Float64,1}:
          0.5278870445904642
          0.5773244171704983
         -0.5241431526595993

        """
function LUPsolve(A, P, b)
    
    #dimension of array
    N = size(A)[1]
    
    #solution vector
    x = Array{Float64, 1}(undef, N)
    
    #permute b same way A was permuted
    b = P*b
    
    #forward substitute L (don't need to use a y, because each row calculation is independent)
    for i = 1:N
        x[i] = b[i]
        for j = 1:(i-1)
            x[i] -= A[i,j]*x[j]
        end
    end

    #backward substitue U
    for i=N:-1:1
        for j=i+1:N
            x[i] -= A[i,j] * x[j]
        end
        x[i] = x[i] / A[i,i]
    end
    return x
end



function conj_grad(A, b, x₀, ϵ::Float64, maxiter::Int)

    x = zeros(length(x₀))
    
    rₚ = 0
    rₙ = b - A*x₀
    p = rₙ
    β = 0

    r_mag = norm(rₙ)
    x_mag = norm(x)

    iter = 0
    while (r_mag/x_mag) > ϵ && iter < maxiter

        p = rₙ + β*p
        α = (rₙ'*rₙ)/(p'*A*p)
        x = x + α * p
        
        rₚ = rₙ
        rₙ = b - A*x
        β = (rₙ'*rₙ)/(rₚ'*rₚ)

        r_mag = norm(rₙ)
        x_mag = norm(x)

        iter += 1
        
    end
    
    return x
end

end # module
