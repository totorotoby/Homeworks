from numpy.random import exponential
import numpy

samples = []
for _ in range(1000000):
    samples.append(exponential(scale = .7))
    

print(numpy.var(samples))
print(numpy.mean(samples))
