function [x_k_plus_1_minus,P_k_plus_1_minus,K_k_plus_1] = EKF_predict(x_k_plus,P_k_plus)
% Takes in the previous state and covariance estimate and returns the
% predicted state,covariance and Gain matrix

% Load Q and R
load cooplocalization_finalproj_KFdata.mat Qtrue Rtrue;

% Perform a state update
x_k_plus_1_minus=EKF_get_next_state(x_k_plus);

% Get F and Omega matrices
[F_tilde,Omega_tilde]=EKF_get_tilde_matrices(x_k_plus);

% Get output matrix
H_tilde=EKF_get_ouput_matrix(x_k_plus_1_minus);

% Propagate Covariance matrix
P_k_plus_1_minus=F_tilde*P_k_plus*F_tilde' + Omega_tilde*Qtrue*Omega_tilde';

% Compute Gain matrix
K_k_plus_1=(P_k_plus_1_minus*H_tilde')/(H_tilde*P_k_plus_1_minus*H_tilde'+Rtrue);

end

