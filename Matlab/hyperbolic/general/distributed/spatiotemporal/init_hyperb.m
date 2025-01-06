
% parameters initialization for the hyperbolic system

% length of the system
L=50;         

% system eigenvalues (characteristic velocities)
%lambda1 = 1;      % always > 0
%lambda2 = 0.02;    % > 0 for collocated inputs

lambda1 = 1;      % always > 0
lambda2 = 0.2;    % > 0 for collocated inputs

% elements of the coupling matrix

 k11 = -0.05;
 k12 = 0.05;
 k21 = 0.05;
 k22 = -0.05;

K = [k11, k12; k21, k22];
LAMBDA = diag([lambda1, lambda2]);

N = 10;