%% Lunar lander example W/observer

% Open the simulink file before starting



close all

%% Define parameters

m=100;
g = 1.62;

rg = 0.5;
j = m*rg^2;
r = 1;

%Main thruster saturation
T_low  = 0;
T_high = 2*m*g;

%Side thruster saturation
S_low  = -0.1*m*g;
S_high =  0.1*m*g;



sc =        [2 2 0.1 0.5 0.5 0.1]'; %initial error scaling
intensity = 0*0.001*[50 50 1 1]';   %measurement noise intensity
bias =      0*[0 0 0 0.1]';           %measurement bias


%% Make system

%States
%x_1: x
%x_2: y
%x_3: pitch

%u_1: side thruster
%u_2: main thruster

A = [zeros(3),eye(3);[0 0 g;0 0 0;0 0 0],zeros(3)];
B = [zeros(3,2);[1/m 0;0 1/m; r/j 0]];

C = [eye(4) zeros(4,2)];
%C = [zeros(4,1) eye(4) zeros(4,1)];

D = zeros(4,2);


sys = ss(A,B,C,D);


rank(ctrb(sys));
rank(obsv(sys));

%% Give cost
M = [eye(3) zeros(3)];

Q = diag([1 1 100]);
R = diag([1 0.001]);
N = zeros(3,2);


%% LQR

syslqr = ss(A,B,M,zeros(3,2));

[K,S,e] = lqry(syslqr,Q,R,N);

syscl = ss(A-B*K,B,C,D);


Kcheck = inv(R)*B'*S; %Check original expression


%% Reference feedforward

H0=evalfr(syscl,0);

P = inv(H0(1:2,1:2));

sysrf = ss(A-B*K,B*P,C,D);


%% Observer

es = eig(A-B*K);

r0 = max(abs(es));

% Radial multiplier & sector
fr = 5;
phi = pi/8;
r   = r0*fr;

spread = -phi:(phi/(2.5)):phi;

poles = -r*exp(1i*spread);


plot(real(es),imag(es),'or',real(poles),imag(poles),'kx');grid on; axis equal
 

L = place(A',C',poles).';

%% Simulate


sim('rocketsim_observer')

%% Plot

tim = nresp.time;

xl  = linresp.signals.values(:,1);
yl  = linresp.signals.values(:,2);
thl = linresp.signals.values(:,3);

xn  = nresp.signals.values(:,1);
yn  = nresp.signals.values(:,2);
thn = nresp.signals.values(:,3);

u1l = linresp_inp.signals.values(:,1);
u2l = linresp_inp.signals.values(:,2);

u1n = nresp_inp.signals.values(:,1);
u2n = nresp_inp.signals.values(:,2);



rx = resp_ref.signals.values(:,1);
ry = resp_ref.signals.values(:,2);


xm = measured_sigs.signals.values(:,1);
xt = true_sigs.signals.values(:,1);

ym = measured_sigs.signals.values(:,2);
yt = true_sigs.signals.values(:,2);


X =     [ -0.5 0.0 0.5  0.5  0.75  0.80  0.70  0.75  0.5  0.0  0.2 -0.2  0.0  -0.5 -0.75 -0.80 -0.70 -0.75 -0.5 -0.5];
Y = 0.7*[  1.0 1.5 1.0 -1.0 -1.60 -1.60 -1.60 -1.60 -1.0 -1.0 -1.4 -1.4 -1.0  -1.0 -1.60 -1.60 -1.60 -1.60 -1.0  1.0];

nT = size(tim,1);

figure 
ani = plot(X,Y,':r','linesmoothing','on');
hold on
plot([-50 200],[0,0],'k','linewidth',2)
hold on
nani = plot(X,Y,'b','linesmoothing','on');
hold on
traj = plot([0,0],[0,0],'.r','linesmoothing','on');
hold on
ntraj = plot([0,0],[0,0],'.b','linesmoothing','on');
hold off
axis equal
axis([-5 20 -5 15])

for k = 1:nT
    Xth = xl(k) + X*cos(thl(k)) + Y*sin(thl(k));
    Yth = yl(k) + -X*sin(thl(k)) + Y*cos(thl(k));
    set(ani,'Xdata',Xth,'Ydata',Yth)

    str = sprintf('Time: %.2f',tim(k));
    title(str)
    
    Xnth = xn(k) + X*cos(thn(k)) + Y*sin(thn(k));
    Ynth = yn(k) + -X*sin(thn(k)) + Y*cos(thn(k));
    set(nani,'Xdata',Xnth,'Ydata',Ynth)
    

    drawnow
end

    
    set(traj,'Xdata',xl(1:15:k),'Ydata',yl(1:15:k))
    set(ntraj,'Xdata',xn(1:15:k),'Ydata',yn(1:15:k))

figure 


subplot 211
plot(tim,180*thl/pi,'--r',tim,xl,'--k',tim,yl,'--b',tim,180*thn/pi,'r',tim,xn,'k',tim,rx,':k',tim,yn,'b',tim,ry,':b','linesmoothing','on')
axis([0 60 -20 30])
legend('\theta','x','y')
title('Plant vs. estimator')

subplot 212
plot(tim,180*(thn-thl)/pi,'r',tim,(xn-xl),'k',tim,(yn-yl),'b','linesmoothing','on')
legend('\theta_e','x_e','y_e')



figure


subplot 211
plot(tim,u1n,'k','linesmoothing','on')
axis([0 60 -0.2*m*g 0.2*m*g])
legend('u_1')
title('Control inputs')


subplot 212
plot(tim,u2n,'k','linesmoothing','on')
axis([0 60 -2*m*g 2*m*g])
legend('u_2')
xlabel('t')


figure

subplot 211
plot(tim,xt,'k',tim,xm,'--r','linesmoothing','on')
legend('x','x_m')
title('True vs. measured')



subplot 212
plot(tim,yt,'k',tim,ym,'--r','linesmoothing','on')
legend('y','y_m')
xlabel('t')



 