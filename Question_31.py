import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def Fuselage_drag_coeff(fus_width, fus_length, fus_height, wings_width, wings_thickness):
    Ld = fus_length/fus_width

    FF = 1 + 2.2/(Ld*1.5)**1.5 + 3.8/(Ld)**3  #form factor

    Re = 1e6 #Reynolds number of 1 million seems reasonable DOES IT????? WHERE HAVE YOU GOT THIS FROM????

    Cf = (1/(3.46*np.log10(Re)-5.6))**2 #skin friction coefficient

    Wet_area = (fus_length*fus_width + fus_height*fus_length + fus_height*fus_width) *2* 2/3 #2/3 correction for the tail getting narrower

    Fus_area = fus_width*fus_height + wings_width*wings_thickness

    C_D_fus = Cf *FF *Wet_area/Fus_area #fuselage drag

    return C_D_fus, Wet_area

def func(CT1, lambda_i, D, W, omega, V):
    CT2 = 2*lambda_i*np.sqrt( (V/(omega*R))**2 + (V*(D/W)/(omega*R) + lambda_i)**2)
    return CT1-CT2

#-------------------------- STANDARD PARAMETERS -------------------------
Vmax = 170*0.5144 #m/s max

rho = 1.225  # Density at sea level (kg/m^3)
cla = 2*np.pi  # lift curve slope [1/deg]
g = 9.81 # m/s^2 gravitational acceleration
mass = 4536 # kg mass of the helicopter MTOW
W = mass * g # N helicopter weight
diam = 13.41 # m Rotor diameter
R = diam/2 # m Rotor radius 
c = 0.76 # m chord length
rotor_RPM = 324 
rotor_speed = rotor_RPM * 2 * np.pi/60  # rad/s angular rotor speed
N = 2 # Number of blades

# Physical quanities of the helicopter
fus_width = 0.99 #m
fus_length  = 13.77 #m
fus_height = 3.11 #m
wings_width = 3.56-0.99 #m the wing width is the total width of the wings minus the fuselage width
wings_thickness = 0.1 #m an approximation of the wing thickness
tail_distance = 7.901 #m, estimated based on the distance of the main rotor on a picture form vertipedia

Cdf, S = Fuselage_drag_coeff(fus_width, fus_length, fus_height, wings_width, wings_thickness)

# Blade solidity
sigma = N * c / (np.pi * R)  # Solidity ratio

#------------------- QUESTION 3 TRIMMING ---------------------
V = np.arange(0, Vmax, 1) # m/s velocity

a1 = []
theta0 = []
for i in range(len(V)):

    # PARAMETER DEFINITIONS
    mu = V[i]/(rotor_speed*R)
    D = 0.5 * rho * V[i]**2 * S * Cdf # drag force N

    CT1 = np.sqrt(W**2+D**2)/(rho*((rotor_speed*R)**2)*np.pi*R**2)

    lambda_i_guess = 0.05
    #CT1, lambda_i, D, W, omega, V
    lambda_i = fsolve(lambda x: func(CT1, x, D, W, rotor_speed, V[i]), lambda_i_guess)

    A = np.array([
        [(1+(3/2)*mu**2), (-(8/3)*mu)],
        [(-mu), ((1+(2/3)*mu**2))]
        ])

    B = np.array([
        -2*(mu)* lambda_i[0] - 2*(mu**2)*(D/W),
        mu*(D/W) + lambda_i[0] + 4*CT1/(cla*sigma)
         ])

    X = np.linalg.solve(A, B) # solve the system of equations for a_1 and theta_0
    
    a1.append(X[0]*180/np.pi) # convert to degrees
    theta0.append(X[1]*180/np.pi) # convert to degrees

plt.plot(V, a1, label='Cyclic input', color = 'orange')
plt.plot(V, theta0, label='Collective input', color = 'purple')
plt.grid()
plt.legend()
plt.xlabel('Velocity (m/s)')
plt.ylabel(r'Pitch Angle (\degree)')
plt.title('Cyclic and Collective Pitch vs Velocity')
plt.show()
