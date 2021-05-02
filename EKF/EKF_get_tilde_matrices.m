function [F_tilde,Omega_tilde] = EKF_get_tilde_matrices(x_k_plus)
% Takes in the state at time step k and returns the F and Omega matrices

% Declare control inputs
v_g=2;
v_a=12;

% Breakout state vector
[eta_g,nu_g,theta_g,eta_a,nu_a,theta_a]=EKF_breakout_state_vector(x_k_plus);

% Construct F matrix using I+dt*A
F_tilde=[1 0 -0.1*v_g*sin(theta_g) 0 0 0;...
         0 1 0.1*v_g*cos(theta_g) 0 0 0;...
         0 0 1 0 0 0;...
         0 0 0 1 0 0.1*v_a*sin(theta_a);...
         0 0 0 0 1 0.1*v_a*cos(theta_a);...
         0 0 0 0 0 1];

% Construct Omega matrix using dt*Gamma      
Omega_tilde=0.1*ones(6,1);
     
end

