function [x_hist,P_hist] = EKF_store_estimate(x_k_plus,P_k_plus)

% Takes in a corrected state and covariance estimate and returns 
% the state and covariance history

% Declare as global for updates in an arbitrary scope
global x_hist P_hist;

% Append corrected state and covariance estimates
x_hist(6,end+1)=x_k_plus;
P_hist(6,6,end+1)=P_k_plus;
end

