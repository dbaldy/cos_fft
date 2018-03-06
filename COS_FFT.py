import math
import cmath
import matplotlib.pyplot as plt
import numpy as np

#COS-FFT Project

# Some random initial variables for testing purposes. They should be later autocalculated
# a and b are the bounderies of the integration
a = -4
b = 4
# k is the iteration in the sum of the integral, N is max(k)
k = 0
N = 50
#Parameters of the distribution
mu = 0
sigma = 2

# characteristic function of Normal Distribution (note: other distributions should be added at a later point)

def CharNormDistr(u,a,b,mu,sigma):
    CharFun =cmath.exp(complex(0,u*mu)-(1/2)*u**2*sigma**2)
    return CharFun

# algo as in the paper
def ProbFunction(x, N,a,b,mu,sigma):
# note: there might be a way to optimise this algo s.t. it calculates F_K once only.
    Probability = 0
    weight = 0.5
    for k in range(0,N-1):
        F_k = (2/(b-a))*(CharNormDistr((k*math.pi)/(b-a),a,b,mu,sigma)*cmath.exp(complex(0,-((k*a*math.pi)/(b-a))))).real
        Probability = Probability + weight*math.cos(k*math.pi*(x-a)/(b-a))*F_k
        weight=1
    return Probability

x = np.arange(a, b, 0.1)
ProbFunctionv= np.vectorize(ProbFunction)
plt.plot(x, ProbFunctionv(x,N,a,b,mu,sigma))

# dx = (b-a)/N
# xmid = numpy.linspace(a+dx/2,b-dx/2,N)
# for k in range(1, N):
#         numpy.strip
# for k = 1:N:
# strip(k) = f(xmid(k))*dx;
# end
# integral = sum(strip)

# Plotting TBA
