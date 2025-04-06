clear all;
clc;


% ----------------------------------------------------------------------------------------
%INITIAL DATA HELICOPTER
g=9.81;	
cla=2*pi; %Airfoil theory [1/rad]
volh=.0722;	%blade solidity	
lok=9.334;
cds=2.415;
mass=4536;
W = mass*g;
rho=1.225;
diam=13.41;
omega = 324*2*pi/60;
vtip= omega*(diam/2);
iy=(1/12)*mass*(4.1^2 + 16.14^2); % Rectangular box approximation with a = 4.1 and b = 16.41 [kg/m^4]
mast=2.285; %(m) Vertical distance between rotor CG and rotor hub
area=pi/4*diam^2;
tau=.1;		%time constant in dynamiCs inflow!!!
collect0=6*pi/180; % (rad) Collective pitch
longit0=0*pi/180; % (rad) Longitudinal cyclic pitch angle
[Cdf, S] = fuselage_dragS();

% ---------------------------------Parameters------------------------------------
%Initial values;
t0=0; %(sec) Setting initial time 
u0=90*0.51444; %(m/sec) Setting initial helicopter airspeed component along body x-axis 
w0=0; %(m/sec) Setting initial helicopter airspeed component along body z-axis
q0=0; %(rad/sec) Setting initial helicopter pitch rate 
pitch0=0*pi/180; %(rad) Setting initial helicopter pitch angle
x0=0; %(m) Setting initial helicopter longitudinal position 
V0 = sqrt(u0^2 + w0^2);
labi0=lambda_i(V0, omega, diam/2, rho, S, Cdf, W); %Initialization non-dimensional inflow in hover!!!

% Angle of attack 
phi=atan(w0/u0);
alfc=longit0-phi;

% mu and lambda_c
qdiml=q0/omega;
vdiml=sqrt(u0^2+w0^2)/vtip;
mu=vdiml*cos(alfc); %Vcos(alphac)/omega*R
labc=vdiml*sin(alfc); %Vsin(alphac)/omega*R

%a1 Flapping calculi
teller=-16/lok*qdiml+8/3*mu*collect0-2*mu*(labc+labi0);
a1=teller/(1-.5*mu^2);

%Thrust coefficient 
ctelem=cla*volh/4*(2/3*collect0*(1+1.5*mu^2)-(labc+labi0));

%Thrust coefficient from Glauert
alfd=alfc-a1; % alpha_d
ctglau=2*labi0*sqrt((vdiml*cos(alfd))^2+(vdiml*sin(alfd)+labi0)^2);
D = 0.5*rho*(V0^2)*cds*S;

labidot=ctelem; 
thrust=labidot*rho*vtip^2*area;

helling=longit0-a1;
vv=vdiml*vtip; %it is 1/sqrt(u^2+w^2)

%===================================================
% Linearised state space
%===================================================
syms u w q pitch collective longit real
syms collect longit real    % Input variables

% Define your nonlinear ODEs symbolically
udot = -g*sin(pitch)-cds/mass*.5*rho*u*vv+...
    thrust/mass*sin(helling)-q*w; % Your full udot
wdot = g*cos(pitch)-cds/mass*.5*rho*w*vv-...
    thrust/mass*cos(helling)+q*u; % Your wdot equation
qdot =-thrust*mast/iy*sin(helling);
    % Your qdot equation
pitchdot = q;

% State and input vectors
x = [u; w; pitch;q];
u_input = [collect; longit];

% Compute Jacobians
A = jacobian([udot; wdot; pitchdot; qdot], x);
B = jacobian([udot; wdot; pitchdot; qdot], u_input);

% Convert to numeric matrices at trim condition
A_num = double(subs(A, {u, w, pitch, q, collect, longit}, {u0, w0, pitch0, q0, collect0, longit0}));
B_num = double(subs(B, {u, w, pitch, q, collect, longit}, {u0, w0, pitch0, q0, collect0, longit0}));



disp(A_num);
D = [0, 0];
Cq = [0, 0, 0, 1];
sys = ss(A_num, B_num, Cq, D);

input_idx = 2;  
output_idx = 4; % Pitch angle output

% Get transfer function
G = tf(sys(output_idx, input_idx));
C = 1;

rlocus(-G*C);






