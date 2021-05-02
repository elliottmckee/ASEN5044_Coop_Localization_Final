function [dxdt] = EKF_dynamics(t_span, x_0)

% Takes in time span and initial condition and returns dxdt to ODE45 solver

% Initialise control inputs
v_g=2;
phi_g=-pi/18;
v_a=12;
omega_a=pi/25;

% Break out state vector
[eta_g, nu_g, theta_g, eta_a, nu_a, theta_a] = EKF_breakout_state_vector(x_0);

% For the EKF, we assume the state updates are noiseless
eta_g_dot=v_g*cos(theta_g);
nu_g_dot=v_g*sin(theta_g);
theta_g_dot=(v_g/0.5)*tan(phi_g);

eta_a_dot=v_a*cos(theta_a);
nu_a_dot=v_a*sin(theta_a);
theta_a_dot=omega_a;

% Construct state vector

dxdt=[eta_g_dot;...
      nu_g_dot;...
      theta_g_dot;...
      eta_a_dot;...
      nu_a_dot;...
      theta_a_dot];

end

