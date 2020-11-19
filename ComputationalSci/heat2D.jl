#=
Solves:

uₜ = k(Δu) + F(x,y,t)        0 ≤ x ≤ L , 0 ≤ y ≤ L, 0 ≤ t ≤ T


u(0,y,t) = u(L,y,t) = u(x,0,t) = u(x,L,t) = 0
u(x,y,0) = f(x,y)

Using standard centered differences in space and forward euler in time.
=#


using LinearAlgebra
using Plots
using Printf



L = 1.0
T = 2.0
λ = .1





F(x,y,t) = (-2 + 2*π^2
    
