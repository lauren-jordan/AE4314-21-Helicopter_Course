import matplotlib.pyplot as plt
import numpy as np

 #------------------------- QUESTION 1.1 -------------------------
def induced_hover(W, rho, R):
    vi = np.sqrt(W / (2 * rho * np.pi * R**2))
    return vi

def density_ISA(h, rho0, T0):
    # Troposphere 0<=h<11km
    alpha = -0.0065 # K/m
    g = 9.81 # m/s^2
    R = 287 # J/kgK
    
    # Density at altitude h
    rho = rho0*(1 - alpha*h/T0)**(g/(R*alpha) - 1)
    return rho

def power_ideal(W, rho, R):
    Pid = W*(induced_hover(W, rho, R))
    return Pid

if __name__ == "__main__":
    # Plot the graph of induced velocity vs hover altitude (ACT)
    graph_vi = False

    # Plot the graph of ideal power vs hover altitude (ACT)
    graph_pid = False

    # Plot the graph of power vs hover altitude (ACT)
    graph_p = False

    # Print the values of the induced velocity, ideal power and power at sea level and at 2500m (ACT)
    print_hover = False

    #------------------------- STANDARD PARAMETERS -------------------------
    Vcruise = 123*0.5144 #knots maximum speed
    Vmax = 170*0.5144 #knots never exceed speed

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

    #------------------------- QUESTION 1.2 - ACT HOVER -------------------------
    # Calculating the induced velocity and ideal power up to the hover ceiling
    vi_rho = [] #m/s
    pid_rho = [] #MW
    p_rho = [] #MW

    for h in range(0, hover_ceiling):
        rho = density_ISA(h, rho0, T0)
        vi = induced_hover(W, rho, R)
        pid = power_ideal(W, rho, R)
        pid_rho.append(pid/10**6)
        p_rho.append(pid/(FM*10**6))
        vi_rho.append(vi)

    if graph_vi == True:
        plt.plot(range(0, hover_ceiling), vi_rho, color = 'purple')
        plt.grid()
        plt.xlabel('Hover Altitude [m]')
        plt.ylabel('Induced Velocity [m/s]')
        plt.title('Induced Velocity vs Hover Altitude')
        plt.show()
    
    if graph_pid == True:
        plt.plot(range(0, hover_ceiling), pid_rho, color = 'pink')
        plt.grid()
        plt.xlabel('Hover Altitude [m]')
        plt.ylabel('Ideal Power [MW]')
        plt.title('Ideal Power vs Hover Altitude (ACT)')
        plt.show()

    if graph_p == True:
        plt.plot(range(0, hover_ceiling), p_rho, color = 'orange')
        plt.grid()
        plt.xlabel('Hover Altitude [m]')
        plt.ylabel('Power [MW]')
        plt.title('Power vs Hover Altitude (ACT)')
        plt.show()

    if print_hover == True:
        print("The induced velocity at sea level is: ", np.round(induced_hover(W, rho0, R), 2), "m/s")
        print("The induced velocity at 2500m is: ", np.round(induced_hover(W, density_ISA(2500, rho0, T0), R), 2), "m/s")
        print(" ")
        print("The ideal power at sea level is: ", np.round(pid_rho[0], 2), "MW")
        print("The ideal power at 2500m is: ", np.round(pid_rho[-1], 2), "MW")
        print(" ")
        print("The power at sea level is: ", np.round(p_rho[0], 2), "MW")
        print("The power at 2500m is: ", np.round(p_rho[-1], 2), "MW")
        print(" ")
        print("The power loading at sea level is: ", np.round(W/(pid_rho[0]*10**6), 2))
        print("The power loading at 2500m is: ", np.round(W/(pid_rho[-1]*10**6), 2))

    

       