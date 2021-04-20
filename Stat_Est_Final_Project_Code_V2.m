%% ASEN 5033 State Estimation Project- Coop Localization
%{
Team Members:
-Elliott McKee
-Robert S.
-Vishal

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Problem 1
%Geometric Constants
L = 0.5;
% va = 12;
% wa = pi()/25;

%Symbolic Math Variables
syms xg eg tg xa ea ta vg pg va wa delta_t y v x f

% State Vector
x = [xg;      %[m] AGV East Pos (Xi)
       eg;      %[m] AGV East Pos (eta)
       tg;       %[RAD] AGV Orientation
       xa;       %[m] UAV East position (Xi)
       ea;       %[m] UAV East Pos (eta)
       ta];       %[RAD] UAV 

f = [vg*cos(tg),    vg*sin(tg),vg/L*tan(pg),        va*cos(ta),va*sin(ta),      wa];

A = jacobian(f,x);



%%  Problem 2
% At initial point in the nominal trajectory
a= subs(A,[vg, va, tg, ta], [2, 12, pi/2, -pi/2]);


dt = 0.1; % discrete time

% Variable F matrix
F = expm(A.*delta_t); 

% Initial F Matrix
FF = expm(a.*dt);


%% Check observability 

h = [ atan((ea - eg)/(xa -xg)) - tg ; sqrt((xg - xa)^2 + (eg - ea)^2) ; atan((eg - ea)/(xg -xa)) - ta ; xa ; ea];
u = [vg;pg;va;wa];
H = jacobian(h,x);


xi = [10,0,pi/2,-60,0,-pi/2];
% Initial H matrix at initial nominal trajectory point
HH = subs(H,x,xi);

% Initialize Observability
O = zeros( 5*6,6);

% Caclulate Observability
for idx = 1:6
    O((idx-1)*5+1:idx*5,:) = HH*FF^(idx-1);
end

% If we need to do this at each state we can just do for loop. save ranks
% and stabilities


% Obs Rank = 6 n = 6 full column rank, no observability issues
Obs = rank(O);

% Find the Eigen Values of F
[V, D] = eig(F);


