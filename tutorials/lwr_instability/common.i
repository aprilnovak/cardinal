pin_pitch = 1.4732e-2         # meters
pin_diameter = 1.1176e-2      # meters
pellet_diameter = 0.9855e-2   # meters
heated_length = 3.6576        # meters
T_inlet = 565                 # Kelvin
P_outlet = 15.51e6            # Pascal

# estimated from the power of BEAVRS, which is 3411 MWth. This is then divided among
# the number of assemblies (193) and the number of fuel pins per assembly (264). The
# result is then multiplied by 4 since we are modeling 4 pins. This gives a nominal
# total power of 267781 Watts.
# [178520.6, 267781, 401671.5]
power = 267781

# https://www.tandfonline.com/doi/full/10.13182/NSE16-3?needAccess=true#d1e860
# This value is parametrically varied to allow the coolant density to either obtain a large
# variation (low flowrate) or a small variation (high flowrate). The BEAVRS flowrate is
# 61.5e6 kg/hr, which we divide by the number of assemblies (193) and the number of pins
# per assembly (17*17), then multiply by 2*2 (number of pins in our 2x2 array).
# Finally, multiply by 1/3600 to convert from hours to seconds. This gives a
# nominal value of ~ 1.2 kg/s
# values: [0.8, 1.2, 1.8] kg/s
mdot = 10.0

fuel_conductivity = 20.0

# number of layers in the initial mesh
n_layers = 50

# number of cells in OpenMC
n_cell_layers = 50
