clear all;
clc;

%===================================================
% Tuner for the loops using linearised state space system
%===================================================

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
R= diam/2;
omega = 324*2*pi/60;
vtip= omega*(diam/2);
iy=(1/12)*mass*(4.1^2 + 16.14^2); % Rectangular box approximation with a = 4.1 and b = 16.41 [kg/m^4]
mast=2.285; %(m) Vertical distance between rotor CG and rotor hub
area=pi/4*diam^2;
collect_0=4.9*pi/180; % (rad) Collective pitch
longit_0=1.32*pi/180; % (rad) Longitudinal cyclic pitch angle
[Cdf, S] = fuselage_dragS();

% ---------------------------------Parameters------------------------------------
%Initial values;
t0=0; %(sec) Setting initial time 
V0=70*0.51444; %(m/sec) Setting initial helicopter airspeed component along body x-axis 
D0 = 0.5*rho*(V0^2)*Cdf*S;
pitch0 = atan(-D0/W);
u0 = V0*cos(pitch0);
w0 = V0*sin(pitch0); %(m/sec) Setting initial helicopter airspeed component along body z-axis
q0 = 0; %(rad/sec) Setting initial helicopter pitch rate 
x0=0; %(m) Setting initial helicopter longitudinal position 
labi0=lambda_i(V0, omega, diam/2, rho, S, Cdf, W); %Initialization non-dimensional inflow in hover!!!

% =================================
% Linearised state space
%==================================

%% Linearize non linear system equations
syms u w q theta_f labi collect0 longit0

V = sqrt(u^2 + w^2);
alpha_c = longit0 - atan(w/u); 
mu = V/(omega*R)*cos(alpha_c);
labc = V*sin(alpha_c)/(omega*R);
%teller(i)=-16/lok*qdiml(i)+8/3*mu(i)*collect(i)-2*mu(i)*(labc(i)+labi(i));
%a1(i)=teller(i)/(1-.5*mu(i)^2);
a1 = (8/3*mu*collect0 - 2*mu*(labc + labi)-16/lok*q/omega)/(1-0.5*mu^2);

% ctglau(i)=2*labi(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*sin(alfd(i))+labi(i))^2);
% ctelem(i)=cla*volh/4*(2/3*collect(i)*(1+1.5*mu(i)^2)-(labc(i)+labi(i)));
ctglau = 2*labi*sqrt((V/(omega*R)*cos(alpha_c-a1))^2+(V/(omega*R)*sin(alpha_c-a1)+labi)^2);
ctel = 1/4*cla*volh*(2/3*collect0*(1+3/2*mu^2)-(labc+labi));

%labidot(i)*rho*vtip^2*area;
T = ctel*rho*(omega*R)^2*pi*R^2;

%D = 0.5*rho*(V_tot(i)^2)*cds*S;
D = Cdf*1/2*rho*V^2*S;

% Derivatives for the matrix 
X = (T/mass)*sin(longit0-a1)-(D/mass)*(u/V)-g*sin(pitch0);
dX = jacobian(X, [u, w, theta_f, q, labi, collect0, longit0]);
dX_x = double(subs(dX, [u, w, theta_f, q, labi, collect0, longit0], [u0, w0, pitch0, q0, labi0, collect_0, longit_0]));

Z = -(T/mass)*cos(longit0-a1)-(D/mass)*w/V+g*cos(pitch0);
dZ = jacobian(Z, [u, w, theta_f, q, labi, collect0, longit0]);
dZ_x = double(subs(dZ, [u, w, theta_f, q, labi, collect0, longit0], [u0, w0, pitch0, q0, labi0, collect_0, longit_0]));

M = -(T/iy)*mast*sin(longit0-a1);
dM = jacobian(M, [u, w, theta_f, q, labi, collect0, longit0]);
dM_x = double(subs(dM, [u, w, theta_f, q, labi, collect0, longit0], [u0, w0, pitch0, q0, labi0, collect_0, longit_0]));

tau = 0.1;
Lambda = 1/tau*(ctel-ctglau)/(omega*R);
dLamb = jacobian(Lambda, [u, w, theta_f, q, labi, collect0, longit0]);
dLamb_new = double(subs(dLamb, [u, w, q, theta_f, labi, collect0, longit0], ...
    [u0, w0, pitch0, q0, labi0, collect_0, longit_0]));

%% Linear system
A = [dX_x(1) dX_x(2) (dX_x(3)-g*cos(pitch0)) (dX_x(4)-w0) dX_x(5);
    dZ_x(1) dZ_x(2) (dZ_x(3)-g*sin(pitch0)) (dZ_x(4)+u0) dZ_x(5); 
    0 0 0 1 0;
    dM_x(1) dM_x(2) dM_x(3) dM_x(4) dM_x(5);
    dLamb_new(1) dLamb_new(2) dLamb_new(3) dLamb_new(4) dLamb_new(5)];

B = [dX_x(6) dX_x(7);
    dZ_x(6) dZ_x(7);
    0 0;
    dM_x(6) dM_x(7);
    dLamb_new(6), dLamb_new(7)];


disp("X_u");
disp(dX_x(1));
Xu = dX_x(1);
disp(" ");
disp("M_u");
disp(dM_x(1));
Mu = dM_x(1);
disp(" ");
disp("M_q");
disp(dM_x(4));
Mq = dM_x(4);

C = eye(5,5);
mdl = ss(A,B,C,0);

% eigs(A);
% disp(eigs(A));

%% Verifying linear model
t_end = 50;

t = linspace(0,t_end, t_end/0.1+1);
u = [ones(1,10)*1*pi/180, zeros(1,length(t)-10); zeros(1,length(t))];
y = lsim(mdl, u, t);

figure;
subplot(2, 3, 1);
plot(t, y(:, 1));
legend('u')

subplot(2, 3, 2);
plot(t, y(:, 2));
legend('w')

subplot(2, 3, 3);
plot(t, y(:, 3)*180/pi);
legend('q [deg/s]')

subplot(2, 3, 4);
plot(t, y(:, 4)*180/pi);
legend('Pitch [deg]')

subplot(2, 3, 5);
plot(t, y(:, 5)*180/pi);
legend('Pitch [deg]')
% 
% %---------------------------------------------------
% syms lamb_i
% lamb_i = vpasolve(lamb_i*(lamb_i - Xu)*(lamb_i - Mq) * 1/Mu +g == 0);
% 
% disp("Eigenvalues:")
% disp(lamb_i(1))
% disp(lamb_i(2))
% disp(lamb_i(3))
% 
% disp("Frequency:")
% wn1 = abs(lamb_i(1));
% wn2 = abs(lamb_i(2));
% wn3 = abs(lamb_i(3));
% disp(wn1)
% disp(wn2)
% disp(wn3)
% 
% disp("Damping coeffs:")
% disp(-real(lamb_i(1))/wn1)
% disp(-real(lamb_i(2))/wn2)
% disp(-real(lamb_i(3))/wn3)


