
rho = 1.225;  % Air density at sea level 
R = 7.32;  % BladeRadius [m]
vtip = 292 * 2 * np.pi / 60 * R;  % Tip velocity [m/s]   
c = 0.53;  % Chord [m]
N = 4;  % Number of blades
solidity = N * c / (np.pi * R);
liftCurve = 2 * np.pi; % Approximation according to thin airfoil theory
CD = 0.3;  % Fuselage Drag Coefficient  
m = 9525; % Mass [kg]
S = 3.9; % Frontal area [m^2]
g = 9.80665; % Acc. due to gravity [m/s]
W = m*g; % Weight [N]
area = pi*R^2; % Area [m^2]

% Speed range 
Vmax = 293; % 180.1*0.51444
V = 