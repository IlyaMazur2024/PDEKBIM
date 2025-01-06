
% transmission line
% Lp = 5;        % line length
% R = 25;        % 25 Om/m
% L = 2.5E-3;    % 2.5 mH/m
% G = 0.3E-6;    % 0.3 uS/m
% C = 5E-9;      % 5 nF/m


% transportation pipeline
Lp = 1000;    % pipeline length [m]
D = 0.2;      % pipeline diameter [m]
A = pi*D^2/4; % pipeline cross sectional area [m^2]
rho = 1000;   % water density [kg/m^3]
lambda = 0.01; % friction coefficient 
c = 1500;     % pressure wave speed [m/s^2]
pi = 5E5;      % input pressure [Pa]
qi = 27.5;     % input mass flow [kg/s]


L = 1/A;       % hydraulic inductivity [1/m^2]
C = A/c^2;     % hydraulic capacitance [s^2]
R = lambda*qi/(rho*D*A^2);  % hydraulic resistance [1/(m^2*s)]
G = 0;
