function [eta_g,nu_g,theta_g,eta_a,nu_a,theta_a] = EKF_breakout_state_vector(x_k_plus_1_minus)
%Simply breaks out state vector into individual state variables

%AGV States
eta_g=x_k_plus_1_minus(1);
nu_g=x_k_plus_1_minus(2);
theta_g=x_k_plus_1_minus(3);

%UAV States
eta_a=x_k_plus_1_minus(4);
nu_a=x_k_plus_1_minus(5);
theta_a=x_k_plus_1_minus(6);
end

