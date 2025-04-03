%STANDARD PARAMETERS
rho = 1.225;% Density at sea level
T0 = 288.15; % K
W = 12564*9.81; % N
R = 8.585; % m
c = 0.59; % m
rotor_speed = 246*2*pi/60; % rad/s
N = 4; %Number of blades
t = 0.016; % m

%INITIAL DATA HELICOTPER
g=9.81;	
cla=5.9; % Sikorsky SC-2110 & SSC-A09
cds=3.90193; % Drag coefficient of the fuselage m^2
mass=12564; % Mass of the helicopter
vtip=rotor_speed*c; %Blade tip speed
diam=2*R;
iy= (1/3)*(1800*c*R*t)*R^2; % MMOI 
mast = 1; %(m) Vertical distance between rotor CG and rotor hub; 
omega=vtip/(diam/2); % Rotor tip speed
area=pi/4*diam^2; % Rotor area
tau=.1;		% Time constant in dynamics inflow
collect(1)=6*pi/180;
longit(1)=0*pi/180;

volh=N*c/pi*R;	%blade solidity	
lok= rho*cla*c*(R^4)/iy; % Lock number gamma = rhp*cla*c*R^4/I

% INITIAL VALUES SIMULATION
t0 = 0;                                      %(sec) Setting initial time 
u0 = 0;                                      %(m/sec) Setting initial helicopter airspeed component along body x-axis 
w0 = 0;                                      %(m/sec) Setting initial helicopter airspeed component along body z-axis
q0 = 0;                                      %(rad/sec) Setting initial helicopter pitch rate 
pitch0 = 0 * pi / 180;                    %(rad) Setting initial helicopter pitch angle
x0 = 0;                                   %(m) Setting initial helicopter longitudinal position 
labi0 = sqrt(mass*g/(area*2*rho))/vtip;   % Initialization non-dimensional inflow

t = [t0];            % Initialization time
u = [u0];            % Initialization helicopter airspeed component along body x-axis
w = [w0];            % Initialization helicopter airspeed component along body z-axis
q = [q0];            % Initialization helicopter pitch rate
pitch = [pitch0];    % Initialization helicopter pitch angle
x = [x0];            % Initialization helicopter position
labi = [labi0];      % Initialization helicopter induced inflow
z = [0];             % Initialization helicopter vertical position

%START INTEGRATION FOR 40 SECONDS
aantal=400;
teind=40;
stap=(teind-t0)/aantal;

for i=1:aantal 

   if t(i)>=0.5 & t(i)<=1 longit(i)=1*pi/180; %longit cyclic input
   else longit(i)=0*pi/180;
   end
   
   %no LAW FOR COLLECTIVE
c(i)=u(i)*sin(pitch(i))-w(i)*cos(pitch(i));
h(i)=-z(i);
collect(i)=collect(1);    

%Defining the differential equations

%defining the nondimensional notations
qdiml(i)=q(i)/omega;
vdiml(i)=sqrt(u(i)^2+w(i)^2)/vtip;
if u(i)==0 	if w(i)>0 	phi(i)=pi/2;
        else phi(i)=-pi/2;end
else
phi(i)=atan(w(i)/u(i));
end
if u(i)<0
phi(i)=phi(i)+pi;
end
alfc(i)=longit(i)-phi(i);

mu(i)=vdiml(i)*cos(alfc(i));
labc(i)=vdiml(i)*sin(alfc(i));

%a1 Flapping calculus
teller(i)=-16/lok*qdiml(i)+8/3*mu(i)*collect(i)-2*mu(i)*(labc(i)+labi(i));
a1(i)=teller(i)/(1-.5*mu(i)^2);

%the thrust coefficient
ctelem(i)=cla*volh/4*(2/3*collect(i)*(1+1.5*mu(i)^2)-(labc(i)+labi(i)));
%Thrust coefficient from Glauert
alfd(i)=alfc(i)-a1(i);
ctglau(i)=2*labi(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*...
sin(alfd(i))+labi(i))^2);

%Equations of motion
labidot(i)=ctelem(i); 
thrust(i)=labidot(i)*rho*vtip^2*area;
helling(i)=longit(i)-a1(i);
vv(i)=vdiml(i)*vtip; 		%it is 1/sqrt(u^2+w^2)

udot(i)=-g*sin(pitch(i))-cds/mass*.5*rho*u(i)*vv(i)+...
thrust(i)/mass*sin(helling(i))-q(i)*w(i);

wdot(i)=g*cos(pitch(i))-cds/mass*.5*rho*w(i)*vv(i)-...
thrust(i)/mass*cos(helling(i))+q(i)*u(i);

qdot(i)=-thrust(i)*mast/iy*sin(helling(i));

pitchdot(i)=q(i);

xdot(i)=u(i)*cos(pitch(i))+w(i)*sin(pitch(i));

zdot(i)=-c(i);

labidot(i)=(ctelem(i)-ctglau(i))/tau;

%Euler integration
u(i+1)=u(i)+stap*udot(i);
w(i+1)=w(i)+stap*wdot(i);
q(i+1)=q(i)+stap*qdot(i);
pitch(i+1)=pitch(i)+stap*pitchdot(i);
x(i+1)=x(i)+stap*xdot(i);
labi(i+1)=labi(i)+stap*labidot(i);
z(i+1)=z(i)+stap*zdot(i);
t(i+1)=t(i)+stap;
end;

plot(t(1:400),longit*180/pi),xlabel('t (s)'),ylabel('longit grd'),grid,pause;

plot(t,u,t,pitch*180/pi),xlabel('t (s)'),ylabel('u(m)'), ylabel('pitch(deg)'),grid,pause;
plot(t,x),xlabel('t (s)'),ylabel('x(m)'),grid,pause;
plot(t,w),xlabel('t (s)'),ylabel('w(m)'),grid,pause;
plot(t,q),xlabel('t (s)'),ylabel('q(m)'),grid,pause; 
plot(t,labi),xlabel('t (s)'),ylabel('labi(m)'),grid,pause;
plot(t,-z),xlabel('t (s)'),ylabel('h(m)'),grid,pause;
