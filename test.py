import math

from OptionPricing import GBMPutPrice, HestonPutPrice, CallPrice, LevyCharFunction, VGPutPrice, CGMYPutPrice
import numpy as np

"""Test file -- Computes the prices of Heston put and call as in the paper"""


print("European option with the Heston Model")
hestonPut = HestonPutPrice(100, 100, 1, 0,
            1.5768, 0.57510, 0.0398, 0.0175, -0.57110, 0, 500)
print("Maturity of 1")
print("Heston Put: %s" % hestonPut)
print("Heston Call: %s" % CallPrice(hestonPut, 100, 100, 0, 10, 0))

hestonPut = HestonPutPrice(100, 100, 10, 0,
            1.5768, 0.57510, 0.0398, 0.0175, -0.57110, 0, 500)
print("Maturity of 10")
print("Heston Put: %s" % hestonPut)
print("Heston Call: %s" % CallPrice(hestonPut, 100, 100, 0, 10, 0))

print(" ")
print(" *** ")
print(" ")

print("European option with the GBM Model")
print("Strike is 80")
strike = 80
gbmPut = GBMPutPrice(100, strike, 0.1, 0, 0.1, 0, 0.25, 500)
print("GBM Put: %s" % gbmPut)
print("GBM Call: %s" % CallPrice(gbmPut, 100, strike, 0, 0.1, 0.1))


print("Strike is 100")
strike = 100
gbmPut = GBMPutPrice(100, strike, 0.1, 0, 0.1, 0, 0.25, 500)
print("GBM Put: %s" % gbmPut)
print("GBM Call: %s" % CallPrice(gbmPut, 100, strike, 0, 0.1, 0.1))

print("Strike is 120")
strike = 120
gbmPut = GBMPutPrice(100, strike, 0.1, 0, 0.1, 0, 0.25, 500)
print("GBM Put: %s" % gbmPut)
print("GBM Call: %s" % CallPrice(gbmPut, 100, strike, 0, 0.1, 0.1))


print(" ")
print(" *** ")
print(" ")

print("European option with the Variance Gamma Model")

print("Maturity is 0.1")
maturity = 0.1
#vgPut = VGPutPrice(100, 90, 0.1, 0, maturity, 0, 0.12, 500)
vgPut= "Y = 0 issue with gamma(-Y)"
print("VG Put: %s" % vgPut)
#print("VG Call: %s" % CallPrice(vgPut, 100, 90, 0, maturity, 0.1))

print("Maturity is 1")
maturity = 1
#vgPut = VGPutPrice(100, 90, 0.1, 0, maturity, 0, 0.12, 500)
vgPut= "Y = 0 issue with gamma(-Y)"
print("VG Put: %s" % vgPut)
#print("VG Call: %s" % CallPrice(vgPut, 100, 90, 0, maturity, 0.1))

print(" ")
print(" *** ")
print(" ")

print("European option with the CGMY model")

print("Y is 0.5")
Y = 0.5
cgmyPut = CGMYPutPrice(100, 100, 0.1, 0, 1, 0, 0, 48,Y) #sigma?
print("CGMY Put: %s" % cgmyPut)
print("CGMY Call: %s" % CallPrice(cgmyPut, 100, 100, 0, 1, 0.1))

print("Y is 1.5")
Y = 1.5
cgmyPut = CGMYPutPrice(100, 100, 0.1, 0, 1, 0, 0, 48,Y) #sigma?
print("CGMY Put: %s" % cgmyPut)
print("CGMY Call: %s" % CallPrice(cgmyPut, 100, 100, 0, 1, 0.1))

print("Y is 0.5")
Y = 1.98
cgmyPut = CGMYPutPrice(100, 100, 0.1, 0, 1, 0, 0, 48,Y) #sigma?
print("CGMY Put: %s" % cgmyPut)
print("CGMY Call: %s" % CallPrice(cgmyPut, 100, 100, 0, 1, 0.1))
