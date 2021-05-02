% Clean slate
clear;
clc;
clf;

% Initial estimates
x_0=[];
P_0=[];

% Variables for EKF filter
xnew=x_0;
Pnew=P_0;

xold=[];
Pold=[];

% Main loop for EKF filter
for k=1:1000
    
    [xold,Pold,K_tilde]=EKF_predict(xnew,Pnew); % Prediction step
    [xnew,Pnew]=EKF_correct(xold,Pold,K_tilde,k); % Correction step
    EKF_store_estimate(xnew,Pnew); % Store the estimate
    
end