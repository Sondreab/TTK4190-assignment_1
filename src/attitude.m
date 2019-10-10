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
N  = 4000;                    % number of samples. Should be adjusted

% model parameters
m = 180;
r = 2;

eps_0 = [0 0 0]';

eta = sqrt(1-eps'*eps);
q_0 = [eta eps']';

I = (m*r^2)*eye(3);       % inertia matrix
I_inv = inv(I);

A = [0 0 0 1/2 0    0 ;
     0 0 0  0  1/2  0 ;
     0 0 0  0   0  1/2;
     0 0 0  0   0   0 ;
     0 0 0  0   0   0 ;
     0 0 0  0   0   0 ];
 
 B = [    0         0          0   ;
          0         0          0   ;
          0         0          0   ;
      1/(m*r^2)     0          0   ;
          0     1/(m*r^2)      0   ;
          0         0     1/(m*r^2)];

k_p = 20;
k_d = 400;

K_d = k_d*eye(3);

%%Problem 1.2

A_2 = [     0           0            0          1/2           0           0      ;
            0           0            0           0           1/2          0      ;
            0           0            0           0            0          1/2     ;
       -k_p/(m*r^2)     0            0      -k_d/(m*r^2)      0           0      ;
            0      -k_p/(m*r^2)      0           0      -k_d/(m*r^2)      0      ;
            0           0      -k_p/(m*r^2)      0            0     -k_d/(m*r^2) ];
        
poles = eig(A_2);

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = -20*deg2rad;

phi_d = 0*deg2rad;            % initial Euler angles
theta_d = 15*cos(0)*deg2rad;
psi_d = 10*sin(0)*deg2rad;

q = euler2q(phi,theta,psi); % transform initial Euler angles to q


w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,14);        % memory allocation

%% FOR-END LOOP
q_tilde_table = zeros(N,4);
for i = 1:N+1,
   t = (i-1)*h;                  % time
   
   %%Problem 1.5
   phi_d = 0*deg2rad;            % initial Euler angles
   theta_d = 15*cos(0.1*t)*deg2rad;
   psi_d = 10*sin(0.05*t)*deg2rad;
   
   q_d = euler2q(phi_d,theta_d,psi_d);
   
   q_tilde = [q_d(1)*q(1) - (-q_d(2:4)')*q(2:4);
              q_d(1)*q(2:4) + q(1)*(-q_d(2:4)) + Smtrx(-q_d(2:4))*q(2:4)];
   q_tilde_table(i,:) = q_tilde;
   
   %%Problem 1.6 
   Eu_d = [phi_d, theta_d, psi_d];
   Eu_d_dot = [0, -1.5*sin(0.1*t), 0.5*cos(0.05*t)]';
   
   T_inv = [1    0           -sin(theta)   ;
            0 cos(phi)  cos(theta)*sin(phi);
            0 -sin(phi) cos(theta)*cos(phi)];
        
   w_d = T_inv*Eu_d_dot*deg2rad;
   
   w_tilde = w - w_d;
     
   tau = -K_d*w_tilde - k_p*q_tilde(2:4);              % control law

   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau'];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
end 

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5);
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
hold off;
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');

figure (4); clf;
hold on;
plot(t, q_tilde_table(:,2), 'b');
plot(t, q_tilde_table(:,3), 'r');
plot(t, q_tilde_table(:,4), 'g');
hold off;
grid on;
legend({'$\epsilon_1$', '$\epsilon_2$', '$\epsilon_3$'}, 'Interpreter', 'Latex');
title('$\tilde{\epsilon}$', 'Interpreter', 'Latex');
xlabel('time [s]'); 
ylabel('input [Rad]');