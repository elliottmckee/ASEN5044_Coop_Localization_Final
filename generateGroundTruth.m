function [x_gt, y_gt] = generateGroundTruth()
%GENERATEGROUNDTRUTH Generates a set of Ground Truth Data, x_gt and y_gt


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
% Load in Given Data
load('cooplocalization_finalproj_KFdata.mat')

%tvec = 0: dt : 100;
%Qtrue = diag([0.001,0.001,0.01,0.001,0.001,0.01]);
%Rtrue = ???????????????????

%ODE Options
opts = odeset('RelTol',1e-11,'AbsTol',1e-13);

%Generate Noise Vectors 
W = (chol(Qtrue)*(randn(length(tvec),6))');

% Generate Nominal with No-Noise
W_nom = [0;0;0;0;0;0];
[t,x_nom] = ode45(@(t,x) nonlin_Dynamics(t, x, uVec, L,W_nom), tvec, x0, opts);
%Reshape
x_nom = x_nom';

% Now for x_gt ground truth
x_gt = zeros(6,1001);
x_gt(:,1) = x_nom(:,1);
for idx = 1:length(tvec)-1
    opts = odeset('RelTol',1e-11,'AbsTol',1e-13);
    [tt,xcheck] = ode45(@(t,x) nonlin_Dynamics(t, x, uVec, L,W(:,idx)), tvec(idx:idx+1), x_gt(:,idx), opts);
    x_gt(:,idx+1) = xcheck(end,:)';
    
end



%Nominal Measurement
y_nom = zeros(5, length(tvec));
y_gt = zeros(5, length(tvec));
for i = 1:length(tvec)
    y_nom(:, i) = h(x_nom(:,i));
    y_gt(:, i) = h(x_gt(:,i)) + chol(Rtrue)*randn(5,1);
end


%% Deal with Angle Rollover
%States
x_gt(3,:) = wrapToPi(x_gt(3,:));
x_gt(6,:) = wrapToPi(x_gt(6,:));

x_nom(3,:) = wrapToPi(x_nom(3,:));
x_nom(6,:) = wrapToPi(x_nom(6,:));

%% Deal with Angle Rollover
y_nom(1,:) = wrapToPi(y_nom(1,:));
y_nom(3,:) = wrapToPi(y_nom(3,:));

y_gt(1,:) = wrapToPi(y_gt(1,:));
y_gt(3,:) = wrapToPi(y_gt(3,:));




end

