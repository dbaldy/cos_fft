import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

def HestonCharFunction(w, u_0, Maturity, eta, lamb, rho, uBarre, mu):
  D = cmath.sqrt((lamb - complex(0, rho * eta * w)) ** 2 + (w ** 2 + complex(0, w)) * eta ** 2)
  G = (lamb - complex(0, rho * eta * w) - D)/(lamb - complex(0, rho * eta * w)  + D)
  hestonCharFunction = cmath.exp(complex(0, w * mu * Maturity) + (u_0)/(eta ** 2) *
                             (1 - np.exp(-D * Maturity))/(1 - G * cmath.exp(-D * Maturity)) * (lamb - complex(0, rho * eta * w) - D)) * cmath.exp((lamb * uBarre)/(eta ** 2) * (Maturity * (lamb - complex(0, rho * eta * w) - D) - 2 * cmath.log((1 - G * cmath.exp(-D * Maturity))/(1 - G))))
  return hestonCharFunction


def Chi(k, c, d, a, b):
    return 1 / (1 + ((k * math.pi) / (b - a)) ** 2) * (math.cos(k * math.pi * (d - a)/(b - a))* np.exp(d) -
            math.cos(k * math.pi * (c - a)/(b - a)) * np.exp(c) + ((k * math.pi) / (b - a)) *
            math.sin(k * math.pi * (d - a)/(b - a)) * np.exp(d) - ((k * math.pi)/(b - a)) *
            math.sin(k * math.pi * (c - a)/(b - a)) * np.exp(c))

def Psi(k, c, d, a, b):
    if k == 0:
        return d - c
    else:
        return (b - a)/(k * math.pi) * (math.sin(k * math.pi * (d - a)/(b - a)) - math.sin(k * math.pi * (c - a)/(b - a)))

def HestonCallPrice(x, strike, current_price, maturity, interest_rate, lamb, eta, uBarre, u_0, rho, mu, N):
    c_1 = mu * maturity + (1 - math.exp(-lamb * maturity)) * (uBarre - u_0)/(2 * lamb) - (1/2) * uBarre * maturity
    # calculating c2
    c_2 = 1/(8 * lamb ** 3) * (eta * maturity * lamb * math.exp(-lamb * maturity) * (u_0 - uBarre) * (8 * lamb * rho - 4 * eta))
    + (lamb * rho * eta * (1 - math.exp(-lamb * maturity)) * (16 * uBarre - 8 * u_0))
    + 2 * uBarre * lamb * maturity * (-4 * lamb * rho * eta + eta ** 2 + 4 * lamb ** 2)
    + eta ** 2 *((uBarre - 2 * u_0) * math.exp(-2 * lamb * maturity) + uBarre * (6 * math.exp(-lamb * maturity) - 7) + 2 * u_0)
    + 8 * lamb ** 2 * (u_0 - uBarre) * (1 - math.exp(-lamb * maturity))
    # establishing the bounds
    a = c_1 - 12 * math.sqrt(abs(c_2))
    b = c_1 + 12 * math.sqrt(abs(c_2))
    tau = maturity
    call_price = 0
    for k in range(0, N - 1):
        u_k = 2 / (b - a) * strike * (Chi(k, 0, b, a, b) - Psi(k, 0, b, a, b))
        call_price += HestonCharFunction((k * math.pi)/(b - a), u_0, maturity, eta, lamb, rho, uBarre, mu) * u_k * cmath.exp(complex(0, k * math.pi)*(x - a)/(b - a))
        if k == 0:
            call_price *= 1 / 2
    return np.real(call_price) * strike * math.exp(-interest_rate * tau)


N = 100
strike = 110
current_price = 100
maturity = 1
interest_rate = 0
mu = 0
x = math.log(current_price/strike)
# HestonCallPrice(x, strike, current_price, maturity, interest_rate, lamb, eta, uBarre, u_0, rho, mu, N):
print(HestonCallPrice(x, strike, current_price, maturity, interest_rate, 1.5768, 0.5751, 0.0398, 0.0175, -0.5711, mu, 100))
