%% 
clear;
close all;

%% Trim conditions
rho = 1.225; %kg/m^3
W = 4536*9.81;
syms lambda_i theta_0 theta_c lambda_c
V0 = 36; % problem is V0 = 0
C_Dfus = 0.1061;
S_fus = 79.3812;
w_mr = 324*2*pi/60;
R_mr = 13.41/2;
c_lalpha = 2*pi;
solidity = .0722;
Ln = 9.334;
m = 4536;
g = 9.81;
Iy = (1/12)*m*(4.1^2 + 16.14^2);
h_mr = 2.285;
D0 = C_Dfus*1/2*rho*V0^2*S_fus;
theta_f0 = atan(-D0/W);
u0 = V0*cos(theta_f0); 
w0 = V0*sin(theta_f0);
T0 = sqrt(W^2 + D0^2);
C_T0 = T0/(rho*(w_mr*R_mr)^2*pi*R_mr^2);
eqn1 = C_T0 == 2*lambda_i*sqrt((V0/(w_mr*R_mr)*cos(D0/W))^2 + (V0/(w_mr*R_mr)*sin(D0/W)+lambda_i)^2);
lambda_i0 = vpasolve(eqn1, lambda_i, [0, 100]);
disp(lambda_i0);
% mu0 = V0/(w_mr*R_mr);
% eqn2 = theta_c == (8/3*mu0*theta_0 - 2*mu0*(lambda_c + lambda_i0))/(1 - 1/2*mu0^2);
% eqn3 = lambda_c == mu0*theta_c + mu0*D0/W;
% eqn4 = C_T0 == c_lalpha*solidity/4*(2/3*theta_0*(1+3/2*mu0^2) - (lambda_c + lambda_i0));
% [theta_00, theta_c0, lambda_c0] = vpasolve([eqn2 eqn3 eqn4], [theta_0 theta_c lambda_c]);
q0 = 0;
theta_00 = 4.9*pi/180;
theta_c0 = 1.32*pi/180;
clear syms
% [theta_00_values, theta_c0_values] = thetaVsV0(linspace(0, 100, 101), C_Dfus, S_fus, solidity, W, w_mr, R_mr, c_lalpha);

%% Linearize non linear system equations

syms u w q theta_f lambda_i theta_0 theta_c
V = sqrt(u^2 + w^2);
alpha_c = theta_c - atan(w/u); % not sure about this
mu = V/(w_mr*R_mr)*cos(alpha_c);
lambda_c = V*sin(alpha_c)/(w_mr*R_mr);
a1 = (8/3*mu*theta_0 - 2*mu*(lambda_c + lambda_i)-16/Ln*q/w_mr)/(1-1/2*mu^2);
C_Tglau = 2*lambda_i*sqrt((V/(w_mr*R_mr)*cos(alpha_c-a1))^2+(V/(w_mr*R_mr)*sin(alpha_c-a1)+lambda_i)^2);
C_TBEM = 1/4*c_lalpha*solidity*(2/3*theta_0*(1+3/2*mu^2)-(lambda_c+lambda_i));
T = C_TBEM*rho*(w_mr*R_mr)^2*pi*R_mr^2;
D = C_Dfus*1/2*rho*V^2*S_fus;

X = T/m*sin(theta_c-a1)-D/m*u/V-g*sin(theta_f0);
dX = jacobian(X, [u, w, q, theta_f, lambda_i, theta_0, theta_c]);
dX_new = double(subs(dX, [u, w, q, theta_f, lambda_i, theta_0, theta_c], ...
    [u0, w0, q0, theta_f0, lambda_i0, theta_00, theta_c0]));

Z = -T/m*cos(theta_c-a1)-D/m*w/V+g*cos(theta_f0);
dZ = jacobian(Z, [u, w, q, theta_f, lambda_i, theta_0, theta_c]);
dZ_new = double(subs(dZ, [u, w, q, theta_f, lambda_i, theta_0, theta_c], ...
    [u0, w0, q0, theta_f0, lambda_i0, theta_00, theta_c0]));

M = -T/Iy*h_mr*sin(theta_c-a1);
dM = jacobian(M, [u, w, q, theta_f, lambda_i, theta_0, theta_c]);
dM_new = double(subs(dM, [u, w, q, theta_f, lambda_i, theta_0, theta_c], ...
    [u0, w0, q0, theta_f0, lambda_i0, theta_00, theta_c0]));

tau = 0.3;
Lambda = 1/tau*(C_TBEM-C_Tglau)/(w_mr*R_mr);
dLambda = jacobian(Lambda, [u, w, q, theta_f, lambda_i, theta_0, theta_c]);
dLambda_new = double(subs(dLambda, [u, w, q, theta_f, lambda_i, theta_0, theta_c], ...
    [u0, w0, q0, theta_f0, lambda_i0, theta_00, theta_c0]));

%% Linear system
% x = [u w q theta_f lambda_i]
% u = [theta_0 theta_c]
A = [dX_new(1) dX_new(2) (dX_new(3)-w0) (dX_new(4)) dX_new(5);
    dZ_new(1) dZ_new(2) (dZ_new(3)+u0) (dZ_new(4)) dZ_new(5);
    dM_new(1) dM_new(2) dM_new(3) dM_new(4) dM_new(5);
    0 0 1 0 0;
    dLambda_new(1) dLambda_new(2) dLambda_new(3) dLambda_new(4) dLambda_new(5)];
B = [dX_new(6) dX_new(7);
    dZ_new(6) dZ_new(7);
    dM_new(6) dM_new(7);
    0 0; 
    dLambda_new(6) dLambda_new(7)];
C = eye(5,5);
mdl = ss(A,B,C,0);

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
legend('pitch [deg/s]')

subplot(2, 3, 4);
plot(t, y(:, 4)*180/pi);
legend('Pitch rate [deg/s]')

subplot(2, 3, 5);
plot(t, y(:, 5));
legend('labi')


% %% Verifying linear model
% t_end = 50;
% t = linspace(0,t_end, t_end/0.1+1);
% u = [ones(1,10)*1*pi/180, zeros(1,length(t)-10); zeros(1,length(t))];
% y = lsim(mdl, u, t);
% 
% t = linspace(0,t_end, t_end/0.1+1);
% u = [ones(1,10)*1*pi/180, zeros(1,length(t)-10); zeros(1,length(t))];
% y = lsim(mdl, u, t);
% 
% figure;
% subplot(2, 3, 1);
% plot(t, y(:, 1));
% legend('u')
% 
% subplot(2, 3, 2);
% plot(t, y(:, 2));
% legend('w')
% 
% subplot(2, 3, 3);
% plot(t, y(:, 3)*180/pi);
% legend('q [deg/s]')
% 
% subplot(2, 3, 4);
% plot(t, y(:, 4)*180/pi);
% legend('Pitch [deg]')
% 
% subplot(2, 3, 5);
% plot(t, y(:, 5));
% legend('labi')

% %% Verifying linear model
% t_end = 50;
% t = linspace(0,t_end, t_end/0.1+1);
% u = [ones(1,10)*5*pi/180, zeros(1,length(t)-10); zeros(1,length(t))];
% y = lsim(mdl, u, t);
% figure;
% plot(t, y(:, 1));
% legend('u')
% 
% figure;
% plot(t, y(:, 2));
% legend('w')
% 
% figure;
% plot(t, y(:, 3)*180/pi);
% legend('q [deg/s]')
% 
% figure;
% plot(t, y(:, 4)*180/pi);
% legend('\theta_f [deg]')
% 
% figure;
% plot(t, y(:, 5));
% legend('\lambda_i')
