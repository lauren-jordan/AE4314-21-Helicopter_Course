
import numpy as np
import matplotlib.pyplot as plt
from Question_12_Forward import vi_BEM

# STANDARD PARAMETERS
rho = 1.225  # Density at sea level (kg/m^3)
T0 = 288.15  # Temperature (K)
W = 12564 * 9.81  # Weight (N)
R = 8.585  # Rotor radius (m)
c = 0.59  # Chord length (m)
rotor_speed = 246 * 2 * np.pi / 60  # Rotor speed (rad/s)
N = 4  # Number of blades
t = 0.016  # Blade thickness (m)

g = 9.81  # Gravitational acceleration (m/s^2)
cla = 5.9  # Lift curve slope (Sikorsky SC-2110 & SSC-A09)
cds = 3.90193  # Drag coefficient of the fuselage (m^2)
mass = 12564  # Mass of the helicopter (kg)
vtip = rotor_speed * c  # Blade tip speed (m/s)
diam = 2 * R  # Rotor diameter (m)

# Moment of Inertia (Iy)
iy = (1/3) * (1800 * c * R * t) * R**2  # Moment of Inertia (kg*m^2)

mast = 1  # Vertical distance between rotor CG and rotor hub (m)
omega = vtip / (diam / 2)  # Rotor tip speed (rad/s)
area = np.pi / 4 * diam**2  # Rotor area (m^2)
tau = 0.1  # Time constant in dynamics inflow

# Initial control inputs
collect = [6 * np.pi / 180] # Collective pitch (radians)
longit = [0 * np.pi / 180]  # Longitudinal cyclic pitch (radians)

# Blade solidity
volh = N * c / (np.pi * R)  # Solidity ratio

# Lock number (Gamma)
lok = rho * cla * c * (R**4) / iy  # Lock number