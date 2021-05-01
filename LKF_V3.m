%% Stat. Est. Final Project Part 2
%% LKF

%Authors:
%Elliott McKee
    %This is me somewhat re-writing the LKF_V2 version authored by Robert S.
%Robert S.
    

%% Housekeeping
clc;clear;close all;


%% Load in Simulated Ground Truth Data
load('cooplocalization_finalproj_KFdata.mat')
load('simulated_groundtruth.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% System Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time Discretization
dt = .1;

%Geometeric
L = .5; %[m] Wheelbase?

%Control Constants
v_g = 2; %[m/s] AGV Ground Velocity
phi_g = -pi/18; %[RAD] AGV Steering Angle

v_a = 12; %[m/s] UAV Air Velocity
omg_a = pi/25; %[RAD/s] UAV Angular Velocity

%Constant Control Vector
uVec = [ v_g;           %[m/s] AGV Linear Velocity
              phi_g;        %[RAD] AGV Steering Angle
              v_a;           %[m/s] UAV Linear Velocity
              omg_a];     %[RAD/s] UAV Angular Rate
          

%% Nominal Trajectory
%Initial Conditions
theta_g0 = pi/2; %[RAD] AGV Initial Orientation
theta_a0 = -pi/2; %[RAD] UAV Initial Orientation

xi_g0 = 10; %[m] AGV Initial East Position
eta_g0 = 0;%[m] AGV Initial North Pos

xi_a0 = -60; %[m] AGV Initial East Position
eta_a0 = 0;%[m] AGV Initial North Pos

% Inital State Vector
x0 = [xi_g0;      %[m] AGV East Pos (Xi)
        eta_g0;      %[m] AGV East Pos (eta)
        theta_g0;       %[RAD] AGV Orientation
        xi_a0;       %[m] UAV East position (Xi)
        eta_a0;       %[m] UAV East Pos (eta)
        theta_a0];       %[RAD] UAV 

%% Nominal Trajectory Functions
theta_g = @(t)  theta_g0 + ( v_g * tan(phi_g) / L ) * t;
theta_a = @(t)  theta_a0 + ( pi / 25 ) * t;
    
%% Required Definitions

% A tilde Matrix, as a function of time
A_tilde = @(t) [ 0   0   -v_g*sin(theta_g(t))    0   0   0;
                        0   0    v_g*cos(theta_g(t))    0   0   0;
                        0   0   0                               0   0   0;
                        0   0   0                               0   0   -v_a*sin(theta_a(t));
                        0   0   0                               0   0   v_a*cos(theta_a(t));
                        0   0   0                               0   0   0]  ;

% B tilde Matrix, as a function of time
B_tilde = @(t) [ cos(theta_g(t))   0                                 0                           0;
                        sin(theta_g(t))    0                                 0                           0;  
                        tan(phi_g)/L     v_g / (L*cos(phi_g)^2)   0                           0 ;                              
                        0                      0                                cos(theta_a(t))        0 ;                      
                        0                       0                               sin(theta_a(t))        0 ;                       
                        0                        0                              0                           1 ]   ;   

               
 % Ftilde Matrix
F_tilde = @(t, dT)  eye(6) + dT * A_tilde(t);

% Gtilde Matrix
G_tilde = @(t, dT)   dT * B_tilde(t);

%H Tilde Matrix
H_tilde = @(X) [ (X(5) - X(2))/((X(5) - X(2))^2 + (X(1)-X(4))^2) , (X(1) - X(4))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , -1 , (X(2) - X(5))/((X(5) - X(2))^2 + (X(4) - X(1))^2 ) , (X(4) - X(1))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , 0 ; ...
                 (X(1) - X(4))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)) , (X(2) - X(5))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)), 0 , (X(4) - X(1))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)) , (X(5) - X(2))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)),0;...
                 (X(5) - X(2))/((X(5) - X(2))^2 + (X(1)-X(4))^2) , (X(1) - X(4))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , 0 , (X(2) - X(5))/((X(5) - X(2))^2 + (X(4) - X(1))^2 ) , (X(4) - X(1))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , -1 ; ...
                                       0                         ,                        0                          , 0 ,                           1                        ,                             0                     ,  0; ...
                                       0                         ,                        0                          , 0 ,                           0                        ,                             1                     ,  0];
                                     

                                   
                                   
                                   
                                   
                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                 
%% LINEARIZED KALMAN IMPLEMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize Vectors
%Total States and Measurements
x_hatP = zeros( 6 , length(x_nom)); %Estimated xHat plus
y_hatP = zeros( 5 , length(y_nom)); %Estimated yHat plus

%Perturbation States and Measurements
dx_hatP = zeros( 6 , length(x_nom)); %Perturbation delta_x plus
dy_hatP = zeros( 5 , length(y_nom)); %Perturbation delta_y plus
%du_tildeP = zeros( 6, length(x_nom) ); %ZERO
 
%Three dimensional matrix containing the Covariance Matrix at each time step on Each PAGE
P_vecP = zeros(6,6, length(tvec));
 
 

 %% Initial Guess
 %Initializing X hat plus as difference between nominal and sim ZERO
 dx_hatP(:,1) = zeros(6,1);
 
 %Initializing P to be somewhat Large
 P_vecP(:,:,1) = 1000* eye(6);
 
 
 %% Kalman Filter Implementation
 %ASSUMING GAMMA = EYE(6)
 
 for ii = 1:length(tvec)-1
     
    % Time Update Step 
    dx_tildeMin = F_tilde(tvec(ii), dt) * dx_hatP(:,ii);
    P_min         =  F_tilde(tvec(ii), dt) *  P_vecP(:,:,ii) *  F_tilde(tvec(ii), dt)' + Qtrue;
    
    % Measurement Update/Correction
    K_kp1 = P_min * H_tilde(x_nom(:,ii+1))' * inv( H_tilde(x_nom(:,ii+1)) * P_min * H_tilde(x_nom(:,ii+1))'  + Rtrue);
    del_y =  y_gt(:, ii+1) - y_nom(:, ii+1);
    P_vecP(:,:,ii+1) = (eye(6) - K_kp1 * H_tilde(x_nom(:,ii+1))) * P_min;
    dx_hatP(:,ii+1) = dx_tildeMin + K_kp1 * (del_y - H_tilde(x_nom(:,ii+1)) * dx_tildeMin);
     
 end
 
 
%% Assemble Total Estimated state
x_hatP = x_nom + dx_hatP;

%% Pull 2-Sigma's from P diagonals
sigVec = zeros(6, length(tvec));

for ii = 1:length(tvec)
    %At each timestep, pull diagonal components of P, take sqrt, multiply by two
    sigVec(:,ii) = 2* sqrt(diag(squeeze(P_vecP(:,:,ii)))) ;
end


%% Calculate Estimation Error
errorVec = x_gt - x_hatP;
 
 
%% Plotting State Estimation Errors
figure()
sgtitle('LKF - State Estimation Errors')
hold on;
for k = 1:6
    ax(k) = subplot(6,1,k);
end

subplot(ax(1));
hold on;
plot(tvec, errorVec(1,:));
plot(tvec, +sigVec(1,:),'-- r');
plot(tvec, -sigVec(1,:),'-- r');
xlabel('time (s)');
ylabel('$e_{\xi_g}$','Interpreter','latex');
legend('Error', '2Sig Bounds')
axis([0 100 -20 20])

subplot(ax(2));
hold on;
plot(tvec, errorVec(2,:));
plot(tvec, +sigVec(2,:),'-- r');
plot(tvec, -sigVec(2,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_g}$','Interpreter','latex');
axis([0 100 -20 20])

subplot(ax(3));hold on;
plot(tvec, errorVec(3,:));
plot(tvec, +sigVec(3,:),'-- r');
plot(tvec, -sigVec(3,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_g}$','Interpreter','latex');
axis([0 100 -2 2])

subplot(ax(4));hold on;
plot(tvec, errorVec(4,:));
plot(tvec, +sigVec(4,:),'-- r');
plot(tvec, -sigVec(4,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\xi_a}$','Interpreter','latex');
axis([0 100 -20 20])

subplot(ax(5));hold on;
% plot(tvec,x_out(5,:));
plot(tvec, errorVec(5,:));
plot(tvec, +sigVec(5,:),'-- r');
plot(tvec, -sigVec(5,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_a}$','Interpreter','latex');
axis([0 100 -5 5])

subplot(ax(6));hold on;
% plot(tvec,x_out(6,:));
plot(tvec, errorVec(6,:));
plot(tvec, +sigVec(6,:),'-- r');
plot(tvec, -sigVec(6,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_a}$','Interpreter','latex');
 axis([0 100 -2 2])


%% Plotting States
figure()
sgtitle('LKF - States')
hold on;
for k = 1:6
    ax(k) = subplot(6,1,k);
end

subplot(ax(1));
hold on;
plot(tvec, x_hatP(1,:));
plot(tvec, x_hatP(1,:) + sigVec(1,:), 'r--');
plot(tvec, x_hatP(1,:) - sigVec(1,:), 'r--');
plot(tvec, x_gt(1,:));
%plot(tvec, +sigVec(1,:),'-- r');
%plot(tvec, -sigVec(1,:),'-- r');
xlabel('time (s)');
ylabel('$e_{\xi_g}$','Interpreter','latex');
legend('Estimated', '2Sig Bounds', '', 'Ground Truth')
axis([0 100 -20 20])

subplot(ax(2));
hold on;
plot(tvec, x_hatP(2,:));
plot(tvec, x_hatP(2,:) + sigVec(2,:), 'r--');
plot(tvec, x_hatP(2,:) - sigVec(2,:), 'r--');
plot(tvec, x_gt(2,:));
%plot(tvec, +sigVec(2,:),'-- r');
%plot(tvec, -sigVec(2,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_g}$','Interpreter','latex');
axis([0 100 -10 10])

subplot(ax(3));hold on;
plot(tvec, x_hatP(3,:));
plot(tvec, x_hatP(3,:) + sigVec(3,:), 'r--');
plot(tvec, x_hatP(3,:) - sigVec(3,:), 'r--');
plot(tvec, x_gt(3,:));
%plot(tvec, +sigVec(3,:),'-- r');
%plot(tvec, -sigVec(3,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_g}$','Interpreter','latex');
axis([0 100 -5 5])

subplot(ax(4));hold on;
plot(tvec, x_hatP(4,:));
plot(tvec, x_hatP(4,:) + sigVec(4,:), 'r--');
plot(tvec, x_hatP(4,:) - sigVec(4,:), 'r--');
plot(tvec, x_gt(4,:));
%plot(tvec, +sigVec(4,:),'-- r');
%plot(tvec, -sigVec(4,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\xi_a}$','Interpreter','latex');
axis([0 100 -110 110])

subplot(ax(5));hold on;
plot(tvec, x_hatP(5,:));
plot(tvec, x_hatP(5,:) + sigVec(5,:), 'r--');
plot(tvec, x_hatP(5,:) - sigVec(5,:), 'r--');
plot(tvec, x_gt(5,:));
%plot(tvec, +sigVec(5,:),'-- r');
%plot(tvec, -sigVec(5,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_a}$','Interpreter','latex');
axis([0 100 -110 110])

subplot(ax(6));hold on;
plot(tvec, x_hatP(6,:));
plot(tvec, x_hatP(6,:) + sigVec(6,:), 'r--');
plot(tvec, x_hatP(6,:) - sigVec(6,:), 'r--');
plot(tvec, x_gt(6,:));
%plot(tvec, +sigVec(6,:),'-- r');
%plot(tvec, -sigVec(6,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_a}$','Interpreter','latex');
axis([0 100 -4 4])



%% Plotting Sigma Bounds
figure()
sgtitle('Covariances with Time')
hold on;

plot(tvec, sigVec(1,:),'-- r');
plot(tvec, -sigVec(1,:),'-- r');
plot(tvec, sigVec(2,:),'-- b');
plot(tvec, -sigVec(2,:),'-- b');
plot(tvec, sigVec(3,:),'-- g');
plot(tvec,-sigVec(3,:),'-- g');
plot(tvec, sigVec(4,:),'-- c');
plot(tvec, -sigVec(3,:),'-- c');
plot(tvec, sigVec(5,:),'-- k');
plot(tvec,-sigVec(5,:),'-- k');
plot(tvec, sigVec(6,:),'-- m');
plot(tvec,-sigVec(6,:),'-- m');
 
 
 
 
 
 
 
 
 