function [x_k_plus_1_plus,P_k_plus_1_plus] = EKF_correct(x_k_plus_1_minus, P_k_plus_1_minus,K_k_plus_1,k)
% Takes in the predicted estimate, Gain matrix and returns the corrected
% state and covariance estimate

% Load output data
load cooplocalization_finalproj_KFdata.mat ydata;

% Get output at time step k+1
y_k_plus_1=ydata(:,k+1);

% Predict output at time step k+1
y_k_plus_1_minus=EKF_get_next_output(x_k_plus_1_minus);

% Get the H matrix at time step k+1
H_tilde=EKF_get_ouput_matrix(x_k_plus_1_minus);

% State correction equation
x_k_plus_1_plus=x_k_plus_1_minus + K_k_plus_1*(y_k_plus_1-y_k_plus_1_minus);

% Covariance correction equation
P_k_plus_1_plus=(eye(6)-K_k_plus_1*H_tilde)*P_k_plus_1_minus;

end

