clear; 

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
tau=.1;		%time constant in dynamics inflow!!!
collect(1)=6.326*pi/180; % (rad) Collective pitch
longit(1)=1.867*pi/180; % (rad) Longitudinal cyclic pitch angle
[Cdf, S] = fuselage_dragS();

% ---------------------------------Parameters------------------------------------
%Trim settings;
t0=0; %(sec) Setting initial time 
q0=0; %(rad/sec) Setting initial helicopter pitch rate 
V0 = 90*0.5144;
x0=0; %(m) Setting initial helicopter longitudinal position 
D0 = 0.5*rho*(V0^2)*Cdf*S;
pitch0 = atan(-D0/W);
u0=V0*cos(pitch0); %(m/sec) Setting initial helicopter airspeed component along body x-axis 
w0 = V0*sin(pitch0);
labi0= lambda_i(V0, omega, diam/2, rho, S, Cdf, W); %Initialization non-dimensional inflow in hover!!!
c0 = u0*sin(pitch0)-w0*cos(pitch0);
pitchdes0 = pitch0;

t(1)=t0;
u(1)=u0;
w(1)=w0;
q(1)=q0;
pitch(1)=pitch0;
x(1)=x0;
labi(1)=labi0;
z(1)=0;
V_tot(1) = V0;
c(1) = c0;
pitchdes(1) = pitchdes0;

% INTEGRATION 
aantal=4800;
teind=480;
stap=(teind-t0)/aantal;
int_errorp(1) = 0;
int_errorc(1) =0;
cdes(1) = 0;
hdes(1) = 0;

%  -------------------------------Start of Simulation------------------------------------
for i = 1:aantal
    if t(i) >= 0.5 && t(i) <= 1
        longit(i) = 1 * pi / 180 +longit(1);
        collect(i) = 0 * pi / 180 + collect(1);
        int_errorp(i) = 0;
        int_errorc(i) =0;

    elseif t(i) >= 360 && t(i) <=360.5
        longit(1) = longit(360);
        collect(1) = collect(360);
        q(1) = q(360);
        hdes(1) = z(360);
        cdes(1) = c(360);
        int_errorc(1) = int_errorc(360);
        int_errorp(1) = int_errorp(360); 
        V_tot(1) = 110*0.5144;
        D = 0.5*rho*(V_tot(1)^2)*Cdf*S;
        pitchdes(1) = atan(-D/W);
        pitch(1) = pitch(360);
    % 
    elseif t(i) >= 240 && t(i) <=240.5
        longit(1) = longit(240);
        collect(1) = collect(240);
        q(1) = q(240);
        hdes(1) = z(240);
        cdes(1) = c(240);
        int_errorc(1) = int_errorc(240); 
        int_errorp(1) = int_errorp(240); 
        V_tot(1) = 90*0.5144;
        D = 0.5*rho*(V_tot(1)^2)*Cdf*S;
        pitchdes(1) = atan(-D/W);
        pitch(1) = pitch(240);
    % % 
    elseif t(i) >= 120 && t(i) <= 120.5
        longit(1) = longit(120);
        collect(1) = collect(120);
        q(1) = q(120);
        hdes(1) = z(120);
        cdes(1) = c(120);
        int_errorc(1) = int_errorc(120); 
        int_errorp(1) = int_errorp(120); 
        V_tot(1) = 70*0.5144;
        D = 0.5*rho*(V_tot(1)^2)*Cdf*S;
        pitchdes(1) = atan(-D/W);
        pitch(1) = pitch(120);

    else
        longit(i) = longit(1);
        collect(i) = collect(1);
        int_errorp(i) = 0;
        int_errorc(i) =0;
    end


   %===================================================
   % Controller always on
   %===================================================
   if t(i)>=2

       % Compute the error term
       error = pitchdes(1)-pitch(i);
       errorc = cdes-c(i);
       
       int_errorp(i) = int_errorp(i-1) + error * stap;
       int_errorc(i) = int_errorc(i-1) + errorc * stap;

       K1 = -0.1; % Proportional term 
       K2 = 0.75;  % Derivative term 
       K3 = -0.02;
       longit(i) = K1*(pitchdes(1)-pitch(i)) + K2*q(i) + K3*int_errorp(i); %PD 

       K4 =  0.03;
       K5 = 0.01;
       K6 = -0.75;
       cdes = K6*(hdes - (z(i)));
       collect(i)= collect(1) + K4*(cdes - c(i)) + K5*int_errorc(i);

   end   
 
    c(i)=u(i)*sin(pitch(i))-w(i)*cos(pitch(i));
    h(i)=-z(i);

   %===================================================
   % Parameters
   %===================================================
    qdiml(i)=q(i)/omega;
    vdiml(i)=sqrt(u(i)^2+w(i)^2)/vtip;
    
    % If u velocity component is 0 
    % phi is the inflow angle!!!
    if u(i)==0 	
        if w(i)>0 	
            phi(i)=pi/2;
        else 
            phi(i)=-pi/2;
        end
    else
        phi(i)=atan(w(i)/u(i));
    end
    
    % Accounting for backwards velocity
    if u(i)<0
        phi(i)=phi(i)+pi;
    end 
    
    % Angle of attack 
    alfc(i)=longit(i)-phi(i);
    
    % mu and lambda_c
    mu(i)=vdiml(i)*cos(alfc(i)); 
    labc(i)=vdiml(i)*sin(alfc(i)); 
    
    %a1 Flapping calculi
    teller(i)=-16/lok*qdiml(i)+8/3*mu(i)*collect(i)-2*mu(i)*(labc(i)+labi(i));
    a1(i)=teller(i)/(1-.5*mu(i)^2);
    
    %Thrust coefficient 
    ctelem(i)=cla*volh/4*(2/3*collect(i)*(1+1.5*mu(i)^2)-(labc(i)+labi(i)));

    %Thrust coefficient from Glauert
    alfd(i)=alfc(i)-a1(i); % alpha_d
    ctglau(i)=2*labi(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*sin(alfd(i))+labi(i))^2);
    
   %===================================================
   % Equations of Motion
   %===================================================
    labidot(i)=ctelem(i); 
    thrust(i)=labidot(i)*rho*vtip^2*area;
    helling(i)=longit(i)-a1(i);
    vv(i)=vdiml(i)*vtip; %it is 1/sqrt(u^2+w^2)

    % udot
    udot(i)=-g*sin(pitch(i))-((Cdf)/mass)*.5*S*rho*u(i)*vv(i)+...
    thrust(i)/mass*sin(helling(i))-q(i)*w(i);

    % wdot
    wdot(i)=g*cos(pitch(i))-(Cdf/mass)*.5*S*rho*w(i)*vv(i)-...
    (thrust(i)/mass)*cos(helling(i))+q(i)*u(i);
    
    %qdot 
    qdot(i)=-thrust(i)*(mast/iy)*sin(helling(i));
    
    % theta_f
    pitchdot(i)=q(i);
    
    % Change in longitudinal position
    xdot(i)=u(i)*cos(pitch(i))+w(i)*sin(pitch(i));
    
    % Change in altitude
    zdot(i)=-c(i);

    % Change in lambda_i
    labidot(i)=(ctelem(i)-ctglau(i))/tau;
    labi(i+1)=labi(i) +stap*labidot(i);
    u(i+1)=u(i) + stap*udot(i);
    w(i+1)=w(i) + stap*wdot(i);

    corrdot(i)=cdes-c(i);

    c(i+1) = c(i) + stap*corrdot(i);
    q(i+1)= q(i) + stap*qdot(i);
    pitch(i+1)= pitch(i) + stap*pitchdot(i);
    x(i+1)=x(i) + stap*xdot(i);
    z(i+1)=z(i) + stap*zdot(i);
    t(i+1)=t(i) + stap;
    V_tot(i+1) = sqrt(u(i+1)^2 + w(i+1)^2);
end;

figure;
subplot(2, 2, 1);
plot(t,pitch*180/pi, 'Color',  [0.4, 0, 0.6]),xlabel('Time (s)'),ylabel('Pitch (deg)'), grid;
subplot(2, 2, 2);
plot(t,q*180/pi,'Color',  [0.4, 0, 0.6]) ,xlabel('Time (s)'),ylabel('Pitch rate q (deg/s)'),grid;
subplot(2, 2, 3);
plot(t,z,'Color',  [0.4, 0, 0.6]),xlabel('Time (s)'),ylabel('Altitude (m)'),grid;
subplot(2, 2, 4);
plot(t,c,'Color',  [0.4, 0, 0.6]),xlabel('Time (s)'),ylabel('Vertical Velocity c (m/s)'), grid, pause;

figure;
subplot(1, 2, 1);
plot(t(1:aantal),collect*180/pi,'Color',  [0.4, 0, 0.6]),xlabel('Time (s)'),ylabel('Collective Input (deg)'),grid;
subplot(1, 2, 2);
plot(t(1:aantal),longit*180/pi,'Color',  [0.4, 0, 0.6]),xlabel('Time (s)'),ylabel('Cyclic input (deg)'),grid,pause;
% 
figure;
subplot(1, 3, 1);
plot(t, u,'Color',  [0.4, 0, 0.6]),xlabel('Time (s)'),ylabel('u (m/s)'),grid;
subplot(1, 3, 2);
plot(t, w,'Color',  [0.4, 0, 0.6]),xlabel('Time (s)'),ylabel('w (m/s)'),grid;
subplot(1, 3, 3);
plot(t, V_tot,'Color',  [0.4, 0, 0.6]),xlabel('Time (s)'),ylabel('Total airspeed V (m/s)'),grid, pause;


