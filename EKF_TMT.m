%% Stat. Est. Final Project Part 2
%% LKF

%Authors:
%Elliott McKee
%This is me somewhat re-writing the LKF_V2 version authored by Robert S.
%Utilizes functions created by Vishal
%Robert S.


%% Housekeeping
%clc;clear;%close all;
clear;

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


%Measurement Model
h = @(X)     [ atan2((X(5) - X(2) ),  ( X(4) - X(1))) - X(3);
                    sqrt((X(1) - X(4))^2 + (X(2) - X(5))^2) ; 
                    atan2((X(2) - X(5)),  (X(1) - X(4))) - X(6) ; 
                    X(4) ; 
                    X(5)];


                
%% TMT Setup
%Number of Runs
N = 250;


%Initialize NEES and NIS Full Solution Vectors
eps_x_TOTAL = zeros(N, length(tvec));
eps_y_TOTAL = zeros(N, length(tvec));

for test = 1:N
    
    %% Generate Ground Truth
    
    %Assume P0
    P0 = diag([6.25 6.252 .03 6.25 6.25 .03]) ;
    %Generate GroundTruth
    [x_gt, y_gt] = generateGroundTruthEKF(x0, P0);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EXTENDED KALMAN IMPLEMENTATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Initialize Vectors
    %Total States and Measurements
    x_hatP = zeros( 6 , length(x_nom)); %Estimated xHat plus
    y_hatP = zeros( 5 , length(y_nom)); %Estimated yHat plus
    
    %Three dimensional matrix containing the Covariance Matrix at each time step on Each PAGE
    P_vecP = zeros(6,6, length(tvec));
    S_vec = zeros(5,5, length(tvec));
    
    
    %% Initial Guess
    %Initializing initial xHat plus to be initial nominal
    x_hatP(:,1) = x_nom(:,1);
    
    %Initializing P to be somewhat Large
    P_vecP(:,:,1) = 1000* eye(6);
    P_vecP(:,:,1) = P0;
    
    
    %% EXTENDED Kalman Filter Implementation
    %ASSUMING GAMMA = EYE(6)
    
    
    %% MATRICES TUNING
    
    %Coarse Tuning
    QKal = Qtrue;
    RKal = Rtrue*.8; 
   
    %Matrix Q Finer Tuning
    QKal(1,1) = QKal(1,1);
    QKal(3,3) = QKal(3,3)/1.5;
    QKal(4,4) = QKal(4,4);
    QKal(5,5) = QKal(5,5);
    QKal(5,5) = QKal(5,5);
    
    %Matrix R Tuning
    RKal(1,1) = RKal(1,1)*8;
    RKal(2,2) = RKal(2,2)/2;
    RKal(3,3) = RKal(3,3)*8;

    
    
    
    
    %ODE45 Options
    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
    
    
    for ii = 1:length(tvec)-1
        
        %% Time Update/Prediction Step
        %Zero Noise ODE45 Integration
        [~,xODERAW] = ode45(@(t,x) nonlin_Dynamics(t, x, uVec, L, zeros(6,1)), tvec(ii:ii+1), x_hatP(:,ii), opts);
        x_hat_kp1Min = xODERAW(end,:)';
        % Deal with Angle Rollover
        x_hat_kp1Min(3) = wrapToPi(x_hat_kp1Min(3));
        x_hat_kp1Min(6) = wrapToPi(x_hat_kp1Min(6));
        
        %Predicted Covariance Update
        P_kp1Min         =  F_tilde(tvec(ii), dt) *  P_vecP(:,:,ii) *  F_tilde(tvec(ii), dt)' + QKal;
        
        
        %% Measurement Update/Correction
        %Nonlinear Measurement Evaluation
        yHat_kp1Min = h(x_hat_kp1Min);
        %Angle Wrap
        yHat_kp1Min(1) = wrapToPi(yHat_kp1Min(1));
        yHat_kp1Min(3) = wrapToPi(yHat_kp1Min(3));
        
        %Save yHat Minus to Vector
        y_hatP(:,ii+1) = yHat_kp1Min; 
        
        %Measurement function Jacobian at Timestep
        Htilde_kp1 = H_tilde(x_hat_kp1Min);
        
        %Nonlinear Measurement Innovation
        eytil_kp1 = y_gt(:,ii+1)  -  yHat_kp1Min;
        
        %Angle Wrap
        eytil_kp1(1) = wrapToPi(eytil_kp1(1));
         eytil_kp1(3) = wrapToPi(eytil_kp1(3));
        
        %Kalman Gain
        K_kp1 = P_kp1Min * Htilde_kp1' * inv( Htilde_kp1 * P_kp1Min *  Htilde_kp1'  + RKal);
        
        %Maintain S matrix
        S_vec(:,:,ii+1) = Htilde_kp1 * P_kp1Min *  Htilde_kp1'  + RKal;
        
        %Update State Estimate and Covariance
        P_vecP(:,:,ii+1) = (eye(6) - K_kp1 * Htilde_kp1) * P_kp1Min;
        x_hatP(:,ii+1) = x_hat_kp1Min + K_kp1 * eytil_kp1;
        % Deal with Angle Rollover
        x_hatP(3, ii+1) = wrapToPi(x_hatP(3, ii+1));
        x_hatP(6, ii+1) = wrapToPi(x_hatP(6, ii+1));
        
        
    end
    
    
    %% Pull 2-Sigma's from P diagonals
    sigVec = zeros(6, length(tvec));
    
    for ii = 1:length(tvec)
        %At each timestep, pull diagonal components of P, take sqrt, multiply by two
        sigVec(:,ii) = 2* sqrt(diag(squeeze(P_vecP(:,:,ii)))) ;
    end
    
    
    %% Calculate Estimation Error
    errorVec_x = x_gt - x_hatP;
    errorVec_y = y_gt - y_hatP;
    
    %Wrap errorVector States
    errorVec_x(3,:) =wrapToPi(errorVec_x(3,:));
    errorVec_x(6,:) = wrapToPi(errorVec_x(6,:));
    
    errorVec_y(1,:) =wrapToPi(errorVec_y(1,:));
    errorVec_y(3,:) = wrapToPi(errorVec_y(3,:));
    
    %%  Calculate NEES and NIS at each step
    eps_x = zeros(1,length(tvec));
    eps_y = zeros(1,length(tvec));
    
    for ii = 1:length(tvec)
        eps_x(ii) = errorVec_x(:,ii)' * inv(P_vecP(:,:,ii)) * errorVec_x(:,ii);
        if(ii > 1)
            eps_y(ii) = errorVec_y(:,ii)' * inv(S_vec(:,:,ii)) * errorVec_y(:,ii);
        end
    end
    
    %Add Epsilon Vector from given test to TOTAL epsilon error
    eps_x_TOTAL(test,:) = eps_x;
    eps_y_TOTAL(test,:) = eps_y;
    
end

%% Average Across Test Runs at each Timestep
NEESfinal = mean(eps_x_TOTAL, 1);
NISfinal = mean(eps_y_TOTAL, 1);


%% Determine NEES/NIS Failure Rates
%NEES Failure Rates
r1NEES = chi2inv(.05/2, 6*N)/N;
r2NEES = chi2inv(1-.05/2, 6*N)/N;

%Failure Rate Determination code CREDIT TO MARTIN GRABAU
NEES_fails_low = sum(NEESfinal<r1NEES);
NEES_fails_high = sum(NEESfinal>r2NEES);
NEES_Failure_rate = (NEES_fails_low+NEES_fails_high)/(length(tvec))



%NIS Failure Rates
r1NIS = chi2inv(.05/2, 5*N)/N;
r2NIS = chi2inv(1-.05/2, 5*N)/N;

%Failure Rate Determination code CREDIT TO MARTIN GRABAU
NIS_fails_low = sum(NISfinal<r1NIS);
NIS_fails_high = sum(NISfinal>r2NIS);
NIS_Failure_rate = (NIS_fails_low+NIS_fails_high)/(length(tvec))


%% NEES Plotting

figure()
hold on

plot(tvec, NEESfinal, 'o')
yline(r1NEES, 'r--')
yline(r2NEES,'r--')

xlabel('Time')
ylabel('NEES Statisctic')
title('NEES Testing')

%% NIS Plotting
r1NIS = chi2inv(.05/2, 5*N)/N;
r2NIS = chi2inv(1-.05/2, 5*N)/N;

figure()
hold on

plot(tvec, NISfinal, 'o')
yline(r1NIS, 'r--')
yline(r2NIS,'r--')

xlabel('Time')
ylabel('NIS Statistic')
title('NIS Testing')




%% Plotting State Estimation Errors
figure()
sgtitle('EKF - State Estimation Errors')
hold on;
for k = 1:6
    ax(k) = subplot(6,1,k);
end

subplot(ax(1));
hold on;
plot(tvec, errorVec_x(1,:));
plot(tvec, +sigVec(1,:),'-- r');
plot(tvec, -sigVec(1,:),'-- r');
xlabel('time (s)');
ylabel('$e_{\xi_g}$','Interpreter','latex');
legend('Error', '2Sig Bounds')
axis([0 100 -20 20])

subplot(ax(2));
hold on;
plot(tvec, errorVec_x(2,:));
plot(tvec, +sigVec(2,:),'-- r');
plot(tvec, -sigVec(2,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_g}$','Interpreter','latex');
axis([0 100 -20 20])

subplot(ax(3));hold on;
plot(tvec, errorVec_x(3,:));
plot(tvec, +sigVec(3,:),'-- r');
plot(tvec, -sigVec(3,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_g}$','Interpreter','latex');
axis([0 100 -2 2])

subplot(ax(4));hold on;
plot(tvec, errorVec_x(4,:));
plot(tvec, +sigVec(4,:),'-- r');
plot(tvec, -sigVec(4,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\xi_a}$','Interpreter','latex');
axis([0 100 -20 20])

subplot(ax(5));hold on;
% plot(tvec,x_out(5,:));
plot(tvec, errorVec_x(5,:));
plot(tvec, +sigVec(5,:),'-- r');
plot(tvec, -sigVec(5,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_a}$','Interpreter','latex');
axis([0 100 -5 5])

subplot(ax(6));hold on;
% plot(tvec,x_out(6,:));
plot(tvec, errorVec_x(6,:));
plot(tvec, +sigVec(6,:),'-- r');
plot(tvec, -sigVec(6,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_a}$','Interpreter','latex');
axis([0 100 -2 2])


%% Plotting States
figure()
sgtitle('EKF - States')
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







%% Plotting Measurement Estimation Errors
figure()
sgtitle('EKF - Measurement Errors')
hold on;
for k = 1:5
    ax(k) = subplot(5,1,k);
end


subplot(ax(1));
hold on;
plot(tvec, errorVec_y(1,:));
plot(tvec, +squeeze(S_vec(1,1,:)) ,'-- r');
plot(tvec, -squeeze(S_vec(1,1,:)) ,'-- r');
xlabel('time (s)');
ylabel('$e_{\xi_g}$','Interpreter','latex');
legend('Error', '2Sig Bounds')
%([0 100 -20 20])

subplot(ax(2));
hold on;
plot(tvec, errorVec_y(2,:));
plot(tvec, +squeeze(S_vec(2,2,:)),'-- r');
plot(tvec, -squeeze(S_vec(2,2,:)),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_g}$','Interpreter','latex');
%axis([0 100 -20 20])

subplot(ax(3));hold on;
plot(tvec, errorVec_y(3,:));
plot(tvec, +squeeze(S_vec(3,3,:)),'-- r');
plot(tvec, -squeeze(S_vec(3,3,:)),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_g}$','Interpreter','latex');
%axis([0 100 -2 2])

subplot(ax(4));hold on;
plot(tvec, errorVec_y(4,:));
plot(tvec, +squeeze(S_vec(4,4,:)),'-- r');
plot(tvec, -squeeze(S_vec(4,4,:)),'-- r');
xlabel('time (s)')
ylabel('$e_{\xi_a}$','Interpreter','latex');
%axis([0 100 -20 20])

subplot(ax(5));hold on;
% plot(tvec,x_out(5,:));
plot(tvec, errorVec_y(5,:));
plot(tvec, +squeeze(S_vec(5,5,:)),'-- r');
plot(tvec, -squeeze(S_vec(5,5,:)),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_a}$','Interpreter','latex');
%axis([0 100 -5 5])









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








