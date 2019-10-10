% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = constant
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and Håkon H. Helgesen

%% USER INPUTS
h = 0.1;                     % sample time (s)
N  = 50000;                    % number of samples. Should be adjusted

% No integral effect
kp = 0.01;
kd = 0.08;
ki = 0;

% With integral effect
kp = 0.09;
kd = 0.22;
ki = 0.00008;

T = 20;
K = 0.1;
b = 0.001;
U = 5;

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

% initial states
x = 0;
y = 100;
psi = 0;
r = 0;
z = 0;

X = [x y psi r z];

x_dot = U;
y_dot = 0;
y_dot_C = 0;
psi_dot = 0;
r_dot = 0;
z_dot = 0;

X_dot = [x_dot y_dot psi_dot r_dot z_dot];

y_int = 0;
u = 5;
v = 0;
 
table = zeros(N+1,6);        % memory allocation

%% FOR-END LOOP

for i = 1:N+1,
   t = (i-1)*h;                  % time
   
   delta = -kp*y-kd*y_dot_C-ki*z;
   if delta > deg2rad*20
       delta = deg2rad*20;
   elseif delta < deg2rad*-20
       delta = deg2rad*-20;
   end
   delta;
   
   x_dot = U*cos(deg2rad*psi);
   y_dot = U*sin(deg2rad*psi);
   y_dot_C = U*psi;
   psi_dot = r;
   r_dot = (K*delta+b-r)/T;
   z_dot = y;
   
   x = x + h*x_dot;
   y = y + h*y_dot;
   psi = psi + h*psi_dot;
   r = r + h*r_dot;
   z = z + h*z_dot;
      
   table(i,:) = [t x y psi r delta];  % store data in table

end 

%% PLOT FIGURES
t       = table(:,1);  
x       = table(:,2); 
y       = table(:,3);
psi     = table(:,4);
r       = rad2deg*table(:,5);
tau     = rad2deg*table(:,6);


figure (1); clf;
hold on;
plot(y, x, 'b');
hold off;
grid on;
title('Position');
xlabel('East [m]'); 
ylabel('North [m]');

figure (2); clf;
hold on;
plot(t, r, 'b');
hold off;
grid on;
legend({'$r$'},'Interpreter','latex');
title('Angular velocity');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure (3); clf;
hold on;
plot(t, tau(:,1), 'r');
plot(t, psi, 'b');
hold off;
grid on;
legend({'$\delta$','$\psi$'},'Interpreter','latex');
title('Control input and heading');
xlabel('time [s]'); 
ylabel('angle [deg]');