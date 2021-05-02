function [dxdt] = nonlin_Dynamics(t, x, u, L,W)
%NONLIN_DYNAMICS ODE45 Function Call to Simulate the Robot Localization System Dynamics

%{
Definitions:
State Vector x = [ xi_g
                          eta_g
                          theta_g
                           xi_a
                          eta_a
                          theta_a]
Control Vector u = [ v_g
                                phi_g
                                v_a
                                omg_a]
%}


%% Break out Control Vector
%Control vector in this case is constant, but non-zero
%AGV Constants
v_g     = u(1); %{m/s} Ground Velocity for Nominal Trajectory
phi_g  = u(2); %[RAD] Steering angle for AGV, Nom Trajectory
%UAV Constants
v_a      = u(3); %[m/s] UAV Velocity for Nominal Trajectory
omg_a = u(4); %[RAD/s] UAV angular velocity for Nominal Trajectory




%% Break out State Vector

% Might need to get rid of the Ws. 
% Qtrue = diag([0.001,0.001,0.01,0.001,0.001,0.01]);
% W = diag(random('Normal',0,Qtrue));
xi_g        =  x(1) ;
eta_g      = x(2) ;
theta_g   = x(3) ;
xi_a        =  x(4) ;
eta_a      =  x(5) ;
theta_a   =  x(6) ;


%% Rates of Change of System States
%AGV State Derivatives
d_xi_g        =  v_g * cos(theta_g) + W(1);
d_eta_g      =  v_g * sin(theta_g) + W(2);
d_theta_g   = (v_g * tan(phi_g)) / L  + W(3);

%UAV State Derivatives
d_xi_a        =  v_a * cos(theta_a) + W(4) ;
d_eta_a      =   v_a * sin(theta_a) + W(5);
d_theta_a   =  omg_a + W(6);


%% Assemble into Vector
dxdt = [d_xi_g;
            d_eta_g;
            d_theta_g;
            d_xi_a;
            d_eta_a;
            d_theta_a];


end