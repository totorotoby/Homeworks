using MatrixMarket
using SparseArrays
using LinearAlgebra


n = 8192
sparident = sparse([i for i in 1:n], [i for i in 1:n], ones(n))
b = [ i for i = 1:n]
a = [ i for i = n:-1:1]



mmwrite("ident.mtx", sparident)

io = open("myb.mtx", "w") do io
    println(io, length(b))
    for x in b
        println(io, x)
    end
end



#@show a
#@show b

