import numpy as np
import matplotlib.pyplot as plt
from Question_12_Forward import vi_BEM, Fuselage_drag_coeff
from Question_12_13_Hover import induced_hover

def flap_angle(a0, a1, b1, psi):

    # Flapping angle
    beta = a0 - a1*np.cos(psi) - b1*np.sin(psi)
    dbeta = a1*np.sin(psi) - b1*np.cos(psi)
    ddbeta = a1*np.cos(psi) + b1*np.sin(psi)
  
    return beta, dbeta

def AOA(pitch, lambda_i, dbeta, r, omega, q, p, psi):

    # Angle of attack calculation
    den = r
    num = (lambda_i + (dbeta*r) - (q/omega)*r*np.cos(psi)- (p/omega)*r*np.sin(psi))
    alpha = pitch - num/den
    return alpha 

#-------------------------- QUESTION 2 Graphs -------------------------
flap_graph = True

twisty = False # True for twisty blade, False for straight blade
#-------------------------- QUESTION 2 Flapping -------------------------
# STANDARD PARAMETERS
rho = 1.225  # Density at sea level (kg/m^3)
cla = 2*np.pi  # lift curve slope [1/deg]
g = 9.81 # m/s^2 gravitational acceleration
mass = 4536 # kg mass of the helicopter MTOW
W = mass * g # N helicopter weight
diam = 13.41 # m Rotor diameter
R = diam/2 # m Rotor radius 
c = 0.76 # m chord length
rotor_RPM = 324 
rotor_speed = -rotor_RPM * 2 * np.pi/60  # rad/s angular rotor speed
N = 2 # Number of blades
t =  0.097 # m blade thickness
FM = 0.76
cds = 2.415
k = 1.15

#### TAIL ROTOR PARAMETERS
D_tail = 2.59 # m tail rotor diameter
Rt = D_tail/2 # m tail rotor radius
ct = 0.30 # m tail rotor chord
kt = 1.4
Omega_tail = 1654 # RPM tail rotor
rotor_speed_t = Omega_tail * 2 * np.pi / 60

# Moment of Inertia (Iy)
h_blade_max = 0.097 # maximum thickness of the blade % of chord
h_blade_avg = h_blade_max/2 # average thickness of the blade % of chord
h_blade = h_blade_avg * c # average thickness of the blade in m

density_blade = 1800 # density of the blade material kg/m^3  (CFRP)
m_blade = density_blade * h_blade * c * R # mass of the blade kg
m_blade = m_blade / 4 # mass reduction of the blade due to it not being a solid block

# Physical quanities of the helicopter
fus_width = 0.99 #m
fus_length  = 13.77 #m
fus_height = 3.11 #m
wings_width = 3.56-0.99 #m the wing width is the total width of the wings minus the fuselage width
wings_thickness = 0.1 #m an approximation of the wing thickness
tail_distance = 7.901 #m, estimated based on the distance of the main rotor on a picture form vertipedia

Cdf, S = Fuselage_drag_coeff(fus_width, fus_length, fus_height, wings_width, wings_thickness)

iy = R**2 * m_blade / 3 # blade moment of inertia kg*m^2

# Blade solidity
volh = N * c / (np.pi * R)  # Solidity ratio

# Lock number (Gamma)
lok = rho * cla * c * (R**4) / iy  # Lock number

#-------------------------- Flapping characteristics -------------------------
gamma = lok # Lock number
collect = 6*np.pi / 180 # Collective pitch (radians)
longcyc = 2*np.pi / 180  # Longitudinal cyclic (radians)
latcyc = 1*np.pi / 180  # Lateral cyclic (radians)
V = 0 # m/s
q = 20*np.pi/180 # rad/s
p = 10*np.pi/180 # rad/s

if twisty:
    twist = -0.175 # radians
else: 
    twist = 0

_, vibar, _, _ = vi_BEM(V, rho, W, R, Cdf, S) # Normalised induced velocity
vi = vibar[0]*induced_hover(W, rho, R) # Induced velocity in forward flight
lambda_i = vi/(rotor_speed*R)

# Fourier coefficients 

a0 = (gamma/8)*(collect - (4/3)*(lambda_i) -(4/5)*twist) 

a1 = (-16*q/(gamma*rotor_speed)+ latcyc + p/rotor_speed)

b1 = -16*p/(rotor_speed*gamma) - longcyc - q/rotor_speed

b1 = b1

print("Coning angle:")
print("a0-a1: ", np.round((a0-a1)*180/np.pi,3))
print("a0+a1: ", np.round((a0+a1)*180/np.pi, 3))

#-------------------------- Flapping angle -------------------------

psi = np.linspace(0, 2*np.pi, 360) # radians

beta_f = [] # radians
dbeta_f = [] # radians
for i in np.arange(0, len(psi)): 
    # gamma, collect, mu, lambda_i, lambda_c, p, q, rotor_speed, longcyc, latcyc, psi 
    beta, dbeta = flap_angle(a0, a1, b1, psi[i])
    beta_f.append(beta)
    dbeta_f.append(dbeta)

if flap_graph:
    plt.plot([i*180/np.pi for i in psi], [i*180/np.pi for i in beta_f], color = "black")
    plt.grid()
    plt.xlabel(r"Azimuth angle $\psi [\degree]$")
    plt.ylabel(r"Flapping angle $\beta [\degree]$")
    plt.title(r"Flapping angle $\beta$ vs Azimuth angle $\psi$")
    plt.show()

#-------------------------- AOA -------------------------
r_var = np.linspace(0, R, 50) # m

AOA_f = np.zeros((len(psi), len(r_var))) # Initialize beta_f array

for i in np.arange(0, len(psi)): 
    for r in np.arange(0, len(r_var)):
        if r_var[r] <= 1.3:
            AOA_f[i, r] = np.NaN 
        else:
            #pitch, lambda_c, lambda_i, dbeta, r, omega, q, p, mu, beta, psi
            pitch = collect + longcyc*np.cos(psi[i]) + latcyc*np.sin(psi[i]) - twist*(r_var[r]/R) # radians
            AOA1 = AOA(pitch, lambda_i, dbeta_f[i], r_var[r]/R, rotor_speed, q, p, psi[i])
            
            AOA_f[i, r] = AOA1*180/np.pi 


