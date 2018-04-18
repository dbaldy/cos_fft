import math
from OptionPricing import GBMPutPrice, HestonPutPrice, CallPrice

hestonPut = HestonPutPrice(100, 100, 10, 0,
                      1.5768, 0.57510, 0.0398, 0.0175, -0.57110, 0, 500)
print("Heston Put: %s" % hestonPut)
print("Heston Call: %s" % CallPrice(hestonPut, 100, 100, 0, 10, 0))

gbmPut = GBMPutPrice(current_price=100, strike=120, interest_rate=0.10,
                     dividend=0, maturity=0.10, mu=0, sigma=0.25, N=500)
print("GBM Put: %s" % gbmPut)
print("GBM Call: %s" % CallPrice(gbmPut, 100, 80, 0, 0.1, 0.1))
