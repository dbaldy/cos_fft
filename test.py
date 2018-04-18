import math
from OptionPricing import GBMCallPrice, HestonCallPrice, PutPrice

heston = HestonCallPrice(100, 100, 1, 0,
                      1.5768, 0.5751, 0.0398, 0.0175, -0.5711, 0, 128)
print(GBMCallPrice(100, 120, 0.1, 0, 0.1, 0, 0.25, 2 ** 6))
print(PutPrice(heston, 100, 0, 1, 100, 0))
