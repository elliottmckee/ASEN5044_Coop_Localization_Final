%% Temporary Driver Script for State Estimation FInal Project
% Elliott McKee
% 4/16/2021


%% Housekeeping
clc;clear;close all;


%% Simulation Setup + Definitions
%Initial State Vector
x0 = [  10;         %[m] Initial AGV East Pos
            0;          %[m] Initial AGV North Pos
            pi/2;       %[m] Initial AGV Orientation
            -60;        %[m] Initial UAV East Pos
            0;            %[m] Initial UAV North Pos
            -pi/2];     %[m] Initial UAV Orientation
        
%Perturbation Vector
x_perturb = [   0;
                       1;
                       0;
                       0;
                       0;
                       0.1];
                   
%Total Initial State Vector
x0_tot = x0 + x_perturb;

        
%Constant Control Vector
uVec = [    2;          %[m/s] AGV Linear Velocity
                -pi/18;     %[RAD] AGV Steering Angle
                12;          %[m/s] UAV Linear Velocity
                pi/25];     %[RAD/s] UAV Angular Rate
            
%Geometric Parameters      
L = 0.5;
            
            
            
%% ODE45 Nonlinear Dynamics Integration
tVec = 0:.1:100;

%ODE Options and Call
opts = odeset('RelTol',1e-11,'AbsTol',1e-13);
[t,x] = ode45(@(t,x) nonlin_Dynamics(t, x, uVec, L), tVec, x0_tot, opts) ; 


%% Deal with Angle Rollover
x(:,3) = wrapToPi(x(:,3));
x(:,6) = wrapToPi(x(:,6));

%% Save Output Matrix
save('nonlinear_simdata', 'x')


%% Plotting
figure()
sgtitle('States vs. Time, Full Nonlinear Dynamics Simulation')

%% AGV States
% AGV East Position
subplot(6,1,1)
hold on
plot(tVec, x(:,1))
xlabel('Time [s]')
ylabel('$\xi_g$ [m] ', 'interpreter', 'latex')
axis([0 100 10 20])

% AGV North Position
subplot(6,1,2)
hold on
plot(tVec, x(:,2))
xlabel('Time [s]')
ylabel('$\eta_g$ [m] ', 'interpreter', 'latex')
axis([0 100 -5 5])

% AGV Orientation
subplot(6,1,3)
hold on
plot(tVec, x(:,3))
xlabel('Time [s]')
ylabel('$\theta_g$ [RAD] ', 'interpreter', 'latex')
axis([0 100 -5 5])


%% UAV  States
% UAV East Position
subplot(6,1,4)
hold on
plot(tVec, x(:,4))
xlabel('Time [s]')
ylabel('$\xi_a$ [m] ', 'interpreter', 'latex')
axis([0 100 -200 200])

% UAV North Position
subplot(6,1,5)
hold on
plot(tVec, x(:,5))
xlabel('Time [s]')
ylabel('$\eta_a$ [m] ', 'interpreter', 'latex')
axis([0 100 -200 200])

% UAV Orientation
subplot(6,1,6)
hold on
plot(tVec, x(:,6))
xlabel('Time [s]')
ylabel('$\theta_a$ [RAD] ', 'interpreter', 'latex')
axis([0 100 -5 5])




                
        
        


