import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

def HestonCharFunction(w, u_0, maturity, eta, lamb, rho, uBarre, mu):
    D = cmath.sqrt((lamb - complex(0, rho * eta * w)) ** 2 + (w ** 2 + complex(0, w)) * eta ** 2)
    G = (lamb - complex(0, rho * eta * w) - D)/(lamb - complex(0, rho * eta * w)  + D)
    return cmath.exp(complex(0, w * mu * maturity) + (u_0)/(eta ** 2) *
                             (1 - np.exp(-D * maturity))/(1 - G * cmath.exp(-D * maturity)) *
                             (lamb - complex(0, rho * eta * w) - D)) * cmath.exp((lamb * uBarre)/(eta ** 2) *
                             (maturity * (lamb - complex(0, rho * eta * w) - D) - 2 *
                              cmath.log((1 - G * cmath.exp(-D * maturity))/(1 - G))))

def LevyCharFunction(w, interest_rate, dividend, maturity, sigma, C, Y, M):
     A = cmath.exp(complex(w * (interest_rate - dividend) * maturity) - 1 / 2 * (w * sigma) ** 2 * maturity)
     B = cmath.exp(maturity * C * math.gamma(-Y) * ((M - complex(w)) ** Y - M ** Y + (G + complex(w)) ** Y -
                                                     G ** Y)) if C > 0 else 1
     return A * B

def GetLevyBounds(c_1, c_2, c_4):
    L = 10
    a = c_1 - L * math.sqrt(c_2 + math.sqrt(c_4))
    b = c_1 + L * math.sqrt(c_2 + math.sqrt(c_4))
    return a, b

def Chi(k, c, d, a, b):
    return 1 / (1 + ((k * math.pi) / (b - a)) ** 2) * (math.cos(k * math.pi * (d - a)/(b - a))* np.exp(d) -
            math.cos(k * math.pi * (c - a)/(b - a)) * np.exp(c) + ((k * math.pi) / (b - a)) *
            math.sin(k * math.pi * (d - a)/(b - a)) * np.exp(d) - ((k * math.pi)/(b - a)) *
            math.sin(k * math.pi * (c - a)/(b - a)) * np.exp(c))

def Psi(k, c, d, a, b):
    if k == 0:
        return d - c
    return (b - a)/(k * math.pi) * (math.sin(k * math.pi * (d - a)/(b - a))
                                        - math.sin(k * math.pi * (c - a)/(b - a)))

def VGPutPrice(current_price, strike, interest_rate, dividend, maturity, mu, sigma, N):
    Y = 0
    C = 1 / nu
    G = theta / sigma ** 2 + math.sqrt(theta ** 2 / sigma ** sigma ** 4 + 2 / (nu * sigma ** 2))
    M = - theta / sigma ** 2 + math.sqrt(theta ** 2 / sigma ** sigma ** 4 + 2 / (nu * sigma ** 2))
    c_1 = (mu + theta) * T
    c_2 = (sigma ** 2 + nu * theta ** 2) * T
    c_4 = 3 * (sigma ** 4 + 2 * thetha ** 4 + nu ** 3 + 4 * (sigma * theta * nu) ** 2) * T
    omega = 1 / nu * math.log(1 - theta * nu - sigma ** 2 * nu / 2) # default log has base e

    a, b = GetLevyBounds(c_1, c_2, c_4)
    price = 0
    for k in range(0, N - 1):
        u_k = 2 / (b - a) * (-Chi(k, a, 0, a, b) + Psi(k, a, 0, a, b))
        price += LevyCharFunction((k * math.pi)/(b - a), interest_rate, dividend, maturity, sigma, C, Y, M) * u_k * \
        cmath.exp(complex(0, k * math.pi * (x - a)/(b - a)))
        if k == 0:
            price /= 2
    return np.real(price) * strike * math.exp(-interest_rate * maturity)


def CGMYPutPrice(current_price, strike, interest_rate, dividend, maturity, mu, sigma, N):
    Y = 1.5
    C = 1
    G = 5
    M = 5
    T = 1
    c_1 = mu * T + C * T * math.gamma(1 - Y) * (M ** (Y - 1) - G ** (Y - 1))
    c_2 = sigma ** 2 * T + C * T  * math.gamma(2 - Y) * (M ** (Y - 2) + G ** (Y - 2))
    c_4 = C * T * math.gamma(4 - Y) * (Y ** (Y - 4) + G ** (Y - 4))
    omega = -C * math.gamma(-Y) * ((M - 1) ** Y - M ** Y + (G + 1) ** Y - G ** Y)

    a, b = GetLevyBounds(c_1, c_2, c_4)
    price = 0
    for k in range(0, N - 1):
        u_k = 2 / (b - a) * (-Chi(k, a, 0, a, b) + Psi(k, a, 0, a, b))
        price += LevyCharFunction((k * math.pi)/(b - a), interest_rate, dividend, maturity, sigma, C, Y, M) * u_k * \
        cmath.exp(complex(0, k * math.pi * (x - a)/(b - a)))
        if k == 0:
            price /= 2
    return np.real(price) * strike * math.exp(-interest_rate * maturity)


def GBMPutPrice(current_price, strike, interest_rate, dividend, maturity, mu, sigma, N):
    x = math.log(current_price / strike)
    c_1 = mu * maturity
    c_2 = sigma ** 2 * maturity
    c_4 = 0
    C = 0
    Y = 0
    M = 0
    a, b = GetLevyBounds(c_1, c_2, c_4)

    price = 0
    for k in range(0, N - 1):
        u_k = 2 / (b - a) * (-Chi(k, a, 0, a, b) + Psi(k, a, 0, a, b))
        price += LevyCharFunction((k * math.pi)/(b - a), interest_rate, dividend, maturity, sigma, C, Y, M) * u_k * \
        cmath.exp(complex(0, k * math.pi * (x - a)/(b - a)))
        if k == 0:
            price /= 2
    return np.real(price) * strike * math.exp(-interest_rate * maturity)

def HestonPutPrice(strike, current_price, maturity, interest_rate, lamb,
                    eta, uBarre, u_0, rho, mu, N):
    x = math.log(current_price / strike)
    c_1 = mu * maturity + (1 - math.exp(-lamb * maturity)) * (uBarre - u_0)/(2 * lamb) - (1/2) * uBarre * maturity
    # calculating c2
    c_2 = 1/(8 * lamb ** 3) * (eta * maturity * lamb * math.exp(-lamb * maturity) * (u_0 - uBarre) * \
    (8 * lamb * rho - 4 * eta)) + (lamb * rho * eta * (1 - math.exp(-lamb * maturity)) * (16 * uBarre - 8 * u_0))
    + 2 * uBarre * lamb * maturity * (-4 * lamb * rho * eta + eta ** 2 + 4 * lamb ** 2)
    + eta ** 2 *((uBarre - 2 * u_0) * math.exp(-2 * lamb * maturity) + uBarre * (6 * math.exp(-lamb * maturity) - 7) + 2 * u_0)
    + 8 * lamb ** 2 * (u_0 - uBarre) * (1 - math.exp(-lamb * maturity))
    # establishing the bounds
    a = c_1 - 12 * math.sqrt(abs(c_2))
    b = c_1 + 12 * math.sqrt(abs(c_2))

    price = 0
    for k in range(0, N - 1):
        u_k = 2 / (b - a) * (-Chi(k, a, 0, a, b) + Psi(k, a, 0, a, b))
        price += HestonCharFunction((k * math.pi)/(b - a), u_0, maturity, eta, lamb, rho, uBarre, mu) * u_k * \
        cmath.exp(complex(0, k * math.pi * (x - a)/(b - a)))
        if k == 0:
            price /= 2
    return np.real(price) * strike * math.exp(-interest_rate * maturity)

def CallPrice(price, current_price, strike, dividend, maturity, interest_rate):
    return np.real(price + current_price * cmath.exp(-dividend * maturity) - strike * cmath.exp(-interest_rate * maturity))


# N = 128
# strike = 100
# current_price = 100
# maturity = 1
# interest_rate = 0
# lamb = 1.5768
# eta = 0.5751
# uBarre = 0.0398
# u_0 = 0.0175
# rho = -0.5711
# mu = 0
# print(HestonCallPrice(strike, current_price, maturity, interest_rate,
#                       lamb, eta, uBarre, u_0, rho, mu, N))
