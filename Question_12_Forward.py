import numpy as np 
import matplotlib.pyplot as plt
from Question_12_13_Hover import induced_hover

def Fuselage_drag_coeff(fus_width, fus_length, fus_height, wings_width, wings_thickness):
    Ld = fus_length/fus_width

    FF = 1 + 2.2/(Ld*1.5)**1.5 + 3.8/(Ld)**3  #form factor

    Re = 1e6 #Reynolds number of 1 million seems reasonable DOES IT????? WHERE HAVE YOU GOT THIS FROM????

    Cf = (1/(3.46*np.log10(Re)-5.6))**2 #skin friction coefficient

    Wet_area = (fus_length*fus_width + fus_height*fus_length + fus_height*fus_width) *2* 2/3 #2/3 correction for the tail getting narrower

    Fus_area = fus_width*fus_height + wings_width*wings_thickness

    C_D_fus = Cf *FF *Wet_area/Fus_area #fuselage drag

    return C_D_fus, Wet_area

def vi_BEM(V, rho0, W, R, Cdf, S):
    #STEP 1 
    # Assumption that gamma = 0 due to forward flight 
    Vel = V

    #Cdf = 0.15
    #S = 3.18*4.71 # Max width of frontal area estimated to be a block using the 
                    # 'Solid wall assumptions' 
    Df = 0.5*rho0*(Vel**2)*Cdf*S # N

    # STEP 2
    # Calculating the angle of attack of the disc
    alphad = np.arcsin(Df/W)
    
    # STEP 3
    # Calculating the velocity compnents of the incoming flow
    V_bar = Vel/induced_hover(W, rho0, R)

    vbarsa = V_bar*np.sin(alphad)
    vbarca = V_bar*np.cos(alphad)

    # STEP 4
    # Solving for the induced velocity numerically
    coeffs = [1, 2*V_bar*np.sin(alphad), (vbarca)**2+(vbarsa)**2, 0, -1]
    coeffs_slow = [1, 0, V_bar**2, 0, -1]
    coeffs_90 = [1, V_bar, -1]

    roots = np.roots(coeffs)
    roots_slow = np.roots(coeffs_slow)
    roots_90 = np.roots(coeffs_90)

    vibar = [r.real for r in roots if np.isreal(r) and r.real > 0]
    vibar_slow = [r.real for r in roots_slow if np.isreal(r) and r.real > 0]
    vibar_90 = [r.real for r in roots_90 if np.isreal(r) and r.real > 0]

    return V_bar, vibar, vibar_slow, vibar_90

#------------------------- QUESTION 1.2 BEM HOVER -------------------------
graph = False

#------------------------- STANDARD PARAMETERS -------------------------
Vcruise = 123*0.5144 #m/s cruise
Vmax = 170*0.5144 #m/s max

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
hover_ceiling =3720 # m hover ceiling
FM = 0.76 # figure of merit

# Physical quanities of the helicopter
fus_width = 0.99 #m
fus_length  = 13.77 #m
fus_height = 3.11 #m
wings_width = 3.56-0.99 #m #### the wing width is the total width of the wings minus the fuselage width
wings_thickness = 0.1 #m an approximation of the wing thickness
tail_distance = 7.901 #m, estimated based on the distance of the main rotor on a picture form vertipedia
#------------------------- QUESTION 1.2 - ACT Forward Flight -------------------------
V = np.arange(0, Vmax, 1) 

Cdf, S = Fuselage_drag_coeff(fus_width, fus_length, fus_height, wings_width, wings_thickness)

vi_bar = []
Vbar = []
vi_bar_slow = []
vi_90 = []

for i in np.arange(0, len(V)):
    V_bar, vibar, vibar_slow, vibar_90 = vi_BEM(V[i], rho0, W, R, Cdf, S)
    
    Vbar.append(V_bar)
    vi_bar.append(vibar)
    vi_bar_slow.append(vibar_slow)
    vi_90.append(vibar_90)

if graph == True:
    fig, ax = plt.subplots(1, 1, figsize = (8, 8))

    ax.plot(Vbar, vi_bar, color = 'black', label = 'Induced Velocity')
    ax.plot(Vbar, [1/Vbar[i] for i in range(0, len(Vbar))], color = '#7F00FF', label = r'High-speed flight $\frac{1}{\bar{V}}$', linestyle = '-.')
    ax.plot(Vbar, vi_90, color = '#FF69B4', label = r'Low-speed flight $\frac{1}{\bar{V} + \bar{v_{i}}}$', linestyle = '-.')
    ax.set_xlabel(r"$\bar{V}$ [m/s]")
    ax.set_ylabel(r"$\bar{v_{i}}$ [m/s]")
    ax.set_title("Induced Velocity vs Forward Flight Velocity")
    ax.set_ylim(0, 1.0)
    ax.grid()
    ax.legend()
    plt.show()

