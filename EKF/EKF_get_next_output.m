function [y_k_plus_1_minus] = EKF_get_next_output(x_k_plus_1_minus)

[eta_g,nu_g,theta_g,eta_a,nu_a,theta_a]=EKF_breakout_state_vector(x_k_plus_1_minus);

y1= wrapToPi(atan2((nu_a-nu_g),(eta_a-eta_g)))-theta_g;

y2=sqrt((eta_g-eta_a)^2+(nu_g-nu_a)^2);

y3=wrapToPi(atan2((nu_g-nu_a),(eta_g-eta_a)))-theta_a;

y4=eta_a;

y5=nu_a;

y_k_plus_1_minus=[y1;y2;y3;y4;y5];


end

