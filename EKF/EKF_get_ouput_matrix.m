function [H_tilde_k_plus_1] = EKF_get_ouput_matrix(x_k_plus_1_minus)
% Takes in the predicted state and returns the H_tilde matrix evaluated at
% the state

% Breakout state vector
[eta_g,nu_g,theta_g,eta_a,nu_a,theta_a]=EKF_breakout_state_vector(x_k_plus_1_minus);

rho=((nu_a-nu_g)^2+(eta_a-eta_g)^2);

% Initialize H_tilde matrix term by term
H11=(nu_a-nu_g)/rho;
H12=(eta_g-eta_a)/rho;
H13=-1;
H14=(nu_g-nu_a)/rho;
H15=(eta_a-eta_g)/rho;
H16=0;

H21=(eta_g-eta_a)/sqrt(rho);
H22=(nu_a-nu_g)/sqrt(rho);
H23=0;
H24=(eta_a-eta_g)/sqrt(rho);
H25=(nu_a-nu_g)/sqrt(rho);
H26=0;

H31=(nu_a-nu_g)/rho;
H32=(eta_g-eta_a)/rho;
H33=0;
H34=(nu_g-nu_a)/rho;
H35=(eta_g-eta_a)/rho;
H36=-1;

H41=0;
H42=0;
H43=0;
H44=1;
H45=0;
H46=0;

H51=0;
H52=0;
H53=0;
H54=0;
H55=1;
H56=0;

% Construct H_tilde matrix
H_tilde_k_plus_1=[H11 H12 H13 H14 H15 H16;...
                  H21 H22 H23 H24 H25 H26;...
                  H31 H32 H33 H34 H35 H36;...
                  H41 H42 H43 H44 H45 H46;...
                  H51 H52 H53 H54 H55 H56;...
                  ];






end

