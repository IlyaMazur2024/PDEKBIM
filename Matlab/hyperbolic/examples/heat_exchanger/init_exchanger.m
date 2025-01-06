
% parameter initialization of the PDE model

L=5;         % length of the exchanger

% fluid constant velocities
vt = 1;      % tube-side fluid velocity [m/s], 
             % vt > 0 for both parallel and counter flow
vs = -0.2;    % shell-side fluid velocity [m/s], 
             % vs > 0 for parallel flow and vs < 0 for counter flow

% physical parameters of the fluids and the exchanger              
ct    = 4200; % specific heat of the tube-side fluid [J/(kg*K)] (water)
cs    = 4200; % specific heat of the shell-side fluid [J/(zc kg*K)] (water)
cw    = 500;  % specific heat of the wall material [J/(kg*K)] (steel)
dti	  = 0.09; % inner diameter of the tube [m]
dto	  = 0.10; % outer diameter of the tube [m]
dsi	  = 0.15; % inner diameter of the shell [m]
ht = 6000;    % heat transfer coefficient [W/(m^2*K)]
hs = 6000;    % heat transfer coefficient [W/(m^2*K)]
rhot  = 1000;  % density of the tube-side fluid [kg/m^3] (water)
rhos  = 1000;  % density of the shell-side fluid [kg/m^3] (water)
rhow  = 7800; % density of the wall material [kg/m^3] (steel)

% generalized parameters of the PDE model
k1 = 4*ht/(rhot*ct*dti);
k2 = 4*dti*ht/(rhow*cw*(dto^2-dti^2));
k3 = 4*dto*hs/(rhow*cw*(dto^2-dti^2));
k4 = 4*dto*hs/(rhos*cs*(dsi^2-dto^2));

k = [k1 k2 k3 k4];
v = [vt vs];

  