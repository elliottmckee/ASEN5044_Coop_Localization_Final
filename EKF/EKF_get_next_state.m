function [x_k_plus_1_minus] = EKF_get_next_state(x_k_plus)
% Takes in state at k, time step k and returns the state at k+1 evolved through
% the full NL, noiseless dynamics.

% Use ODE45 to solve non-linear EOMs
[t,x]=ode45(@EKF_dynamics,[0 0.1],x_k_plus);

% Take the state at the latest time step, perform wrapping and return 
x_k_plus_1_minus=EKF_wrap_state(x(end,:));
end

