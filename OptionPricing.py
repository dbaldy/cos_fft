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
  # return hestonCalc
  #exp(1i*omega*mu*T+u0*((1-exp(-D(omega)*T))./(1-G(omega).*exp(-D(omega)*T))).*(lambda-1i*rho*eta*omega-D(omega))/eta^2) ...
#    .* exp((T*(lambda-1i*rho*eta*omega-D(omega))-2*log((1-G(omega).*exp(-D(omega)*T))./(1-G(omega))))*lambda*u_bar/eta^2);



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

def HestonCallPrice(x, strike, current_price, maturity, interest_rate, lamb,
                    eta, uBarre, u_0, rho, mu, N):
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
    
    call_price = 0
    for k in range(0, N - 1):
        u_k = 2 / (b - a) * (Chi(k, 0, b, a, b) - Psi(k, 0, b, a, b))
        call_price += HestonCharFunction((k * math.pi)/(b - a), u_0, maturity, eta, lamb, rho, uBarre, mu) * u_k * \
        cmath.exp(complex(0, k * math.pi * (x - a)/(b - a)))
        if k == 0:
            call_price /= 2
    return np.real(call_price) * strike * math.exp(-interest_rate * maturity)


N = 128
strike = 100
current_price = 100
maturity = 1
interest_rate = 0
lamb = 1.5768
eta = 0.5751
uBarre = 0.0398
u_0 = 0.0175
rho = -0.5711
mu = 0
x = math.log(current_price / strike)
# HestonCallPrice(x, strike, current_price, maturity, interest_rate, lamb, eta, uBarre, u_0, rho, mu, N):
print(HestonCallPrice(x, strike, current_price, maturity, interest_rate,
                      lamb, eta, uBarre, u_0, rho, mu, N))
