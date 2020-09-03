import numpy as np
import math

# Structural engine parameters
# Connecting rod (m)
connRodLength = 183.2e-3
# Crankshaft radius (m)
crankRadius = 43.15e-3
# Piston diameter (m)
pistonDiameter = 70e-3
# Cylinder volume (m³)
cylClearanceVol = 32.64-6
# Clearance volume height (m)
cylclearanceHeight = cylClearanceVol / (math.pi * pistonDiameter/2) ** 2
# Piston area (m²)
pistonArea = math.pi * (pistonDiameter/2) ** 2
# Crank radius to connecting rod length ratio (m)
crankConnRodRatio = crankRadius / connRodLength
# Max piston displacement (m)
maxPistonDisp = 2 * crankRadius
# Total cylinder displaced volume (m³)
cylinderDispVol = cylClearanceVol + maxPistonDisp*pistonArea
# Clearance volume related area at top dead center (m²)
clearanceVolArea = cylClearanceVol/cylclearanceHeight

# Fuel and mixture properties
# Fuel lower heating value (J/kg) [E100 = 24.65e6; E22 = 38.92e6]
lowerFuelHeatVal = 24.65e6
# Fuel mass injected per cyclinder per cycle (kg/cyc)
fuelMass = 40.03e-6
# Max fuel heat release per cylinder per cycle (J/cyc)
maxFuelHeatRelease = lowerFuelHeatVal*fuelMass