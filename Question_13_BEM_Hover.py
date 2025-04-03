import numpy as np 
import matplotlib.pyplot as plt
from Question_12_13_Hover import induced_hover, density_ISA, power_ideal


#------------------------- QUESTION 1.2 BEM HOVER -------------------------

def p_induced_bem(W, R, rho):
    T = W
    k = 1.15
    return k*T*induced_hover(T, rho, R)

def p_profile_bem(N, omega, c, R, cd, rho):
    sigma = N*c/(np.pi*R)
    Pp = (sigma*cd*rho*np.pi*(R**2)*(omega*R)**3)/8
    return Pp

def CL(W, rho, R, omega, N, c):
    sigma = N*c/(np.pi*R)
    T = W
    CL = 6.6*T/(rho*sigma*np.pi*(R**2)*(omega*R)**2)
    return CL

#------------------------- QUESTION 1.3 FORWARD FLIGHT  -------------------------
# Print the BEM values for hover 
print_BEM = False

#------------------------- STANDARD PARAMETERS -------------------------
T0 = 288.15  # Temperature (K)
rho0 = 1.225  # Density at sea level (kg/m^3)
g = 9.81 # m/s^2 gravitational acceleration
mass = 4536 # kg mass of the helicopter MTOW
W = mass * g # N helicopter weight
diam = 13.41 # m Rotor diameter
R = diam/2 # m Rotor radius 
c = 0.76 # m chord length
rotor_RPM = 324 
rotor_speed = rotor_RPM * 2 * np.pi/60  # rad/s angular rotor speed
N = 2 # Number of blades

#------------------------- QUESTION 1.2 BEM HOVER -------------------------
# Calculating the induced power at sea level and at the ceiling
vi0 = induced_hover(W, rho0, R)

# Calculating the lift coefficient estimate
CLbar0 = CL(W, rho0, R, rotor_speed, N, c)

#   Calculating the induced power at sea level and at the ceiling
p_induced0 = p_induced_bem(W, R, rho0)

# Calculating the profile power at sea level and at the ceiling
cd00 = 0.01 # Taken from cl-cd curve 

p_profile0 = p_profile_bem(N, rotor_speed, c, R, cd00, rho0) 

# Calculating the figure of merit at sea level and at the ceiling
FM0 = p_induced0/(p_induced0 + p_profile0)

# Calculating the total power at sea level and at the ceiling
Phov0 = p_induced0 + p_profile0

if __name__ == "__main__":
    if print_BEM == True:
        print("The lift coefficient at sea level is: ", np.round(CLbar0, 2))

        print(" ")

        print("The figure of merit at sea level is: ", np.round(FM0, 3))

        print(" ")
        print("The power at sea level is: ", np.round(Phov0/(FM0*10**6), 2), "MW")

        print(" ")
        print("The power loading at sea level is: ", np.round(W/Phov0, 2), "N/W")





