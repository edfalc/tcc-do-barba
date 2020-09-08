import numpy as np
import math

def get_mapped_value(xInputVec, yInputVec, interpOrder,xValue):
    poly = np.polyfit(xInputVec, yInputVec, interpOrder)
    return np.polyval(poly, xValue)

def piston_disp(crankAngle, crankRadius, crankConnRodRatio):
    return crankRadius * (1 + 1 / crankConnRodRatio - math.cos(crankAngle * math.pi / 180) - 1 / crankConnRodRatio * math.sqrt(1 - crankConnRodRatio ** 2 * (math.sin(crankAngle * math.pi / 180)) ** 2))

def heat_transfer_coeff(cylinderVol, cylinderPress, mixtureTemp, avgPistonSpeed):
    return 130 * cylinderVol ** (-0.06) * (cylinderPress * 1e-5) ** (0.8) * mixtureTemp ** (-0.4) * (avgPistonSpeed + 1.4) ** 0.8

def cylinder_press(gasConst, specHeatConst, cylinderVol, fuelHeatReleaseRateDelta, heaTransfCylWallRate, fuelMassDelta, cylinderPress, cylinderVolDelta):
    return gasConst / (specHeatConst * cylinderVol) * (fuelHeatReleaseRateDelta - heaTransfCylWallRate + fuelMassDelta - (1 + specHeatConst / gasConst) * cylinderPress * cylinderVolDelta)