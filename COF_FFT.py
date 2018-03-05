import math
import cmath
# import pylab
import numpy

#COS-FFT Project

# Some random initial variables for testing purposes. They should be later autocalculated
# a and b are the bounderies of the integration
a = 0
b = 5
# k is the iteration in the sum of the integral, N is max(k)
k = 0
N = 50
#Parameters of the distribution
mu = 0
sigma = 1
x = 1

# characteristic function of Normal Distribution (note: other distributions should be added at a later point)
def CharNormDistr(x,mu,sigma):
    CharFun = cmath.exp(complex(0,x * mu -1)/(2* x**2* sigma))
    return CharFun

# algo as in the paper
def ProbFunction(x, N, func):
# note: there might be a way to optimise this algo s.t. it calculates F_K once only.
    Probability = 0
    for k in range(1,N+1):
        F_k = 2 / (b - a)*(func(x, mu, sigma) *
                           cmath.exp(complex(0,-(k * a * math.pi/(b - a))))).real
        Probability += math.cos(k*math.pi*(x-a)/(b-a))*F_k
    return Probability

print(ProbFunction(0.5, N, CharNormDistr))

# dx = (b-a)/N
# xmid = numpy.linspace(a+dx/2,b-dx/2,N)
# for k in range(1, N):
#         numpy.strip
# for k = 1:N:
# strip(k) = f(xmid(k))*dx;
# end
# integral = sum(strip)

# Plotting TBA
