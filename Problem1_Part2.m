%% ASEN 5044 Question 1 Part 2 Script
% Elliott McKee
% 4/19/2021

%{
NOTES:

Currently just using the hard coded A and B matrices found by hand

%}


%% Housekeeping
clc;clear;%close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1 Question 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% System Description Constants
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


%% Matrix Definition Functions
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


%% Measurement Equation 
% h = @(X)     [ wrapToPi(atan((X(5) - X(2) )/ ( X(4) - X(1)))) - wrapToPi(X(3)) ;
%                     sqrt((X(1) - X(4))^2 + (X(2) - X(5))^2) ; 
%                     wrapToPi(atan((X(2) - X(5))/ (X(1) - X(4)))) - wrapToPi(X(6)) ; 
%                     X(4) ; 
%                     X(5)];

%Measurement Model
h = @(X)     [ wrapToPi(atan2((X(5) - X(2) ),  ( X(4) - X(1)))) - X(3);
                    sqrt((X(1) - X(4))^2 + (X(2) - X(5))^2) ; 
                    wrapToPi(atan2((X(2) - X(5)),  (X(1) - X(4)))) - X(6) ; 
                    X(4) ; 
                    X(5)];


H_tilde = @(X) [ (X(5) - X(2))/((X(5) - X(2))^2 + (X(1)-X(4))^2) , (X(1) - X(4))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , -1 , (X(2) - X(5))/((X(5) - X(2))^2 + (X(4) - X(1))^2 ) , (X(4) - X(1))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , 0 ; ...
                 (X(1) - X(4))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)) , (X(2) - X(5))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)), 0 , (X(4) - X(1))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)) , (X(5) - X(2))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)),0;...
                 (X(5) - X(2))/((X(5) - X(2))^2 + (X(1)-X(4))^2) , (X(1) - X(4))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , 0 , (X(2) - X(5))/((X(5) - X(2))^2 + (X(4) - X(1))^2 ) , (X(4) - X(1))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , -1 ; ...
                                       0                         ,                        0                          , 0 ,                           1                        ,                             0                     ,  0; ...
                                       0                         ,                        0                          , 0 ,                           0                        ,                             1                     ,  0];
             
 

%% Construct Nominal Trajectory and Measurement Vectors
%Time Vector
tVec = 0: dt : 100;

%ODE Options and Call
opts = odeset('RelTol',1e-11,'AbsTol',1e-13);
[t,x_nom] = ode45(@(t,x) nonlin_Dynamics(t, x, uVec, L), tVec, x0, opts);

%Reshape
x_nom = x_nom';


%Nominal Measurement
y_nom = zeros(5, length(tVec));
for i = 1:length(tVec)
    y_nom(:, i) = h(x_nom(:,i));
end

%% Deal with Angle Rollover
y_nom(1,:) = wrapToPi(y_nom(1,:));
y_nom(3,:) = wrapToPi(y_nom(3,:));





%% Propogate State Forward with Time

%Pre-allocate Perturbation State
x_perturb = zeros(6, length(tVec));
%Pre-allocate Full State
x_full = zeros(6, length(tVec));
%Pre-allocate Measurement Perturb Vector
y_perturb = zeros(5, length(tVec));
%Pre-allocate Measurement FULL Vector
y_full = zeros(5, length(tVec));

%Initial Perturbation
x_p0 = [    0;
                1;
                0;
                0;
                0;
                0.1];

%Intial Values
x_perturb(:,1) = x_p0;
x_full(:,1) = x0;
y_perturb(:,1) = H_tilde(x0) * x_perturb(:, i);

for i = 1:length(tVec)-1

    %Step forward in Time- assumes zero perturbation in control state
    x_perturb(:, i+1) = F_tilde(tVec(i), dt) * x_perturb(:, i) ;
    %Reconstruct Total State, using nominal trajectory calculated previously
    x_full(:, i+1) = x_nom(:, i+1) + x_perturb(:, i+1);

    %Calculate Perturbation Measurements
    y_perturb(:, i+1) = H_tilde(x_nom(:, i+1)) * x_perturb(:, i+1);
    %Calcualte Total Measurements
    y_full(:, i+1) = y_nom(:, i+1) + y_perturb(:, i+1);
end



%% Deal with Angle Rollover
%States
x_full(3,:) = wrapToPi(x_full(3,:));
x_full(6,:) = wrapToPi(x_full(6,:));


%% Load Non-linear Simulation Data
nonlineardata = load('nonlinear_simdata.mat')
x_nonlin = nonlineardata.x';
y_nonlin = nonlineardata.y';


%% TOTAL STATE Plotting
figure()
sgtitle('TOTAL SYSTEM STATES vs. Time, Linearized Dynamics Simulation')

%% AGV States
% AGV East Position
subplot(6,1,1)
hold on
plot(tVec, x_full(1,:))
plot(tVec, x_nonlin(1,:), '--')
xlabel('Time [s]')
ylabel('$\xi_g$ [m] ', 'interpreter', 'latex')
legend('Linearized', 'Nonlinear', 'location', 'northeast')
axis([0 100 10 16])

% AGV North Position
subplot(6,1,2)
hold on
plot(tVec, x_full(2,:))
plot(tVec, x_nonlin(2,:), '--')
xlabel('Time [s]')
ylabel('$\eta_g$ [m] ', 'interpreter', 'latex')
axis([0 100 -2 4])

% AGV Orientation
subplot(6,1,3)
hold on
plot(tVec, x_full(3, :))
plot(tVec, x_nonlin(3,:), '--')
xlabel('Time [s]')
ylabel('$\theta_g$ [RAD] ', 'interpreter', 'latex')
axis([0 100 -3 3])


%% UAV  States
% UAV East Position
subplot(6,1,4)
hold on
plot(tVec, x_full(4, :))
plot(tVec, x_nonlin(4,:), '--')
xlabel('Time [s]')
ylabel('$\xi_a$ [m] ', 'interpreter', 'latex')
axis([0 100 -60 150])

% UAV North Position
subplot(6,1,5)
hold on
plot(tVec, x_full(5, :))
plot(tVec, x_nonlin(5,:), '--')
xlabel('Time [s]')
ylabel('$\eta_a$ [m] ', 'interpreter', 'latex')
axis([0 100 -100 100])

% UAV Orientation
subplot(6,1,6)
hold on
plot(tVec, x_full(6, :))
plot(tVec, x_nonlin(6,:), '--')
xlabel('Time [s]')
ylabel('$\theta_a$ [RAD] ', 'interpreter', 'latex')
axis([0 100 -3 3])




%% PERTURBATION STATE Plotting
figure()
sgtitle('Perturbation States vs. Time, Linearized Dynamics Simulation')

%% AGV States
% AGV East Position
subplot(6,1,1)
hold on
plot(tVec, x_perturb(1,:))
xlabel('Time [s]')
ylabel('$\xi_g$ [m] ', 'interpreter', 'latex')
axis([0 100 -1 1])

% AGV North Position
subplot(6,1,2)
hold on
plot(tVec, x_perturb(2,:))
xlabel('Time [s]')
ylabel('$\eta_g$ [m] ', 'interpreter', 'latex')
axis([0 100 0 2])

% AGV Orientation
subplot(6,1,3)
hold on
plot(tVec, x_perturb(3, :))
xlabel('Time [s]')
ylabel('$\theta_g$ [RAD] ', 'interpreter', 'latex')
axis([0 100 -1 1])


%% UAV  States
% UAV East Position
subplot(6,1,4)
hold on
plot(tVec, x_perturb(4, :))
xlabel('Time [s]')
ylabel('$\xi_a$ [m] ', 'interpreter', 'latex')
axis([0 100 -10 10])

% UAV North Position
subplot(6,1,5)
hold on
plot(tVec, x_perturb(5, :))
xlabel('Time [s]')
ylabel('$\eta_a$ [m] ', 'interpreter', 'latex')
axis([0 100 0 20])

% UAV Orientation
subplot(6,1,6)
hold on
plot(tVec, x_perturb(6, :))
xlabel('Time [s]')
ylabel('$\theta_a$ [RAD] ', 'interpreter', 'latex')
axis([0 100 -1 1])




%% MEASUREMENT STATE Plotting
figure()
sgtitle('LINEARIZED MEASUREMENTS vs. Time, Linearized Dynamics Simulation')

%% AGV States
% 1st Measure
subplot(5,1,1)
hold on
plot(tVec, y_full(1,:))
plot(tVec, y_nonlin(:,1), '--')
xlabel('Time [s]')
ylabel('$\gamma_{ag}$ [RAD] ', 'interpreter', 'latex')
legend('Linearized','NonLinear')
%axis([0 100 -1 1])

% 2nd Measure
subplot(5,1,2)
hold on
plot(tVec, y_full(2,:))
plot(tVec, y_nonlin(:,2), '--')
xlabel('Time [s]')
ylabel('$\rho_{ga}$ [m] ', 'interpreter', 'latex')
%axis([0 100 -1 1])

% 3rd Measure
subplot(5,1,3)
hold on
plot(tVec, y_full(3,:))
plot(tVec, y_nonlin(:,3), '--')
xlabel('Time [s]')
ylabel('$\gamma_{ga}$ [RAD] ', 'interpreter', 'latex')
%axis([0 100 -1 1])

% 4th Measure
subplot(5,1,4)
hold on
plot(tVec, y_full(4,:))
plot(tVec, y_nonlin(:,4), '--')
xlabel('Time [s]')
ylabel('$\xi_a$ [m] ', 'interpreter', 'latex')
%axis([0 100 -1 1])

% 5th Measure
subplot(5,1,5)
hold on
plot(tVec, y_full(5,:))
plot(tVec, y_nonlin(:,5), '--')
xlabel('Time [s]')
ylabel('$\eta_a$ [m] ', 'interpreter', 'latex')
%axis([0 100 -1 1])