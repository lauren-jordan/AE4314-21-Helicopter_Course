
clc; clear; close all;

%-------------------------- STANDARD PARAMETERS -------------------------
Vmax = 170 * 0.5144; % m/s max velocity
rho = 1.225;  % Density at sea level (kg/m^3)
cla = 2 * pi;  % Lift curve slope [1/rad]
g = 9.81; % m/s^2 gravitational acceleration
mass = 4536; % kg mass of the helicopter MTOW
W = mass * g; % N helicopter weight
diam = 13.41; % m Rotor diameter
R = diam / 2; % m Rotor radius 
c = 0.76; % m chord length
rotor_RPM = 324;
rotor_speed = rotor_RPM * 2 * pi / 60;  % rad/s angular rotor speed
N = 2; % Number of blades

% Physical quantities of the helicopter
fus_width = 0.99; % m
fus_length = 13.77; % m
fus_height = 3.11; % m
wings_width = 3.56 - 0.99; % m
wings_thickness = 0.1; % m

tail_distance = 7.901; % m, estimated based on a picture from Vertipedia

% Compute fuselage drag coefficient
Ld = fus_length / fus_width;
FF = 1 + 2.2 / (Ld * 1.5)^1.5 + 3.8 / (Ld)^3; % Form factor
Re = 1e6; % Reynolds number
Cf = (1 / (3.46 * log10(Re) - 5.6))^2; % Skin friction coefficient
Wet_area = (fus_length * fus_width + fus_height * fus_length + fus_height * fus_width) * 2 * 2/3;
Fus_area = fus_width * fus_height + wings_width * wings_thickness;
Cdf = Cf * FF * Wet_area / Fus_area;
S = Wet_area;

% Blade solidity
sigma = N * c / (pi * R);

%------------------- QUESTION 3 TRIMMING ---------------------
V = 0:1:Vmax; % m/s velocity
num_V = length(V);

% Initialize arrays
a1 = zeros(1, num_V);
theta0 = zeros(1, num_V);

for i = 1:num_V
    [a1(i), theta0(i)] = compute_pitch(V(i), rotor_speed, R, rho, S, Cdf, W, cla, sigma);
end
% Plot results

figure;
plot(V, a1*180/pi, 'black', 'LineWidth', 1.5); hold on;
plot(V, theta0*180/pi,'Color',  [0.7, 0.6, 0.9], 'LineWidth', 1.5);
grid on;
legend('Cyclic input', 'Collective input');
xlabel('Velocity (m/s)');
ylabel('Pitch Angle (°)');
%title('Cyclic and Collective Pitch vs Velocity');