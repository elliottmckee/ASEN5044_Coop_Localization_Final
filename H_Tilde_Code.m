

%% Project Code
%

% H_Tilde Matrix

% Vars v_g theta_g eta_x xi_x
%        1    2     3       4     5    6   
% x = [xi_g eta_g theta_g xi_a eta_a theta_a]
% t is the column/time step for nominal trajectory


H_tilde = @(X) [ (X(5) - X(2))/((X(5) - X(2))^2 + (X(1)-X(4))^2) , (X(1) - X(4))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , -1 , (X(2) - X(5))/((X(5) - X(2))^2 + (X(4) - X(1))^2 ) , (X(4) - X(1))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , 0 ; ...
                 (X(1) - X(4))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)) , (X(2) - X(5))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)), 0 , (X(4) - X(1))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)) , (X(5) - X(2))/sqrt(((X(5) - X(2))^2 + (X(1)-X(4))^2)),0;...
                 (X(5) - X(2))/((X(5) - X(2))^2 + (X(1)-X(4))^2) , (X(1) - X(4))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , 0 , (X(2) - X(5))/((X(5) - X(2))^2 + (X(4) - X(1))^2 ) , (X(4) - X(1))/((X(5) - X(2))^2 + (X(4) - X(1))^2) , -1 ; ...
                                       0                         ,                        0                          , 0 ,                           1                        ,                             0                     ,  0; ...
                                       0                         ,                        0                          , 0 ,                           0                        ,                             1                     ,  0]
             
             
    



%(eta_a - eta_g)^2
% (X(5) - X(2))^2

% (xi_a - xi_g)^2
%  (X(4) - X(1))^2


%(eta_g - eta_a)
% (X(2) - X(5))

% (xi_g - xi_a)
%  (X(1) - X(4))