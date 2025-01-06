function [gttw,gtsw,gstw,gssw] = calc_exchanger_par_distr_freq_resp(omega,l,v,k)
%
% This function calculates spatially distributed frequency responses 
% for a double-pipe parallel-flow heat exchanger. 
% The frequency responses are calculated from the irrational transfer functions  
% of the heat exchanger resulting from the following PDEs:
%
% dTt/dt + vt*dTt/dl = k1*(Tw-Tt)		       Tt - tube-side fluid temperature
% dTw/dt             = k2*(Tt-Tw)+k3*(Ts-Tw)   Tw - wall temperature
% dTs/dt + vs*dTs/dl = k4*(Tw-Ts)		       Ts - shell-side fluid temperature
%
% Inputs for vt>0 and vs>0 (boundary conditions for parallel-flow):
% Tt(0,t) = Tti(t), Ts(0,t) = Tsi(t) 
%
% Outputs for vt>0 and vs>0 (parallel-flow):
% Tto(t) = Tt(L,t), Tso(t) = Ts(L,t) 
%
% Function inputs: 
% omega - vector of angular frequencies, 
% l - vector of the spatial positions in [0 L],
%     (the first element in l should be equal to 0 and the last one - to L)
% v - vector of the fluid velocities in PDEs, v = [vt vs],
% k - vector of the constant parameters in PDEs, k = [k1 k2 k3 k4].
%
% Function outputs: 
% matrices of spatially distributed frequency responses for the parallel flow:
% gttw(l,i*omega) = Tt(l,i*omega)/Tt(0,i*omega)
% gtsw(l,i*omega) = Tt(l,i*omega)/Ts(0,i*omega)
% gstw(l,i*omega) = Ts(l,i*omega)/Tt(0,i*omega)
% gssw(l,i*omega) = Ts(l,i*omega)/Ts(0,i*omega)
%
% Krzysztof Bartecki, 2020

vt = v(1);  vs = v(2);
k1  = k(1); k2  = k(2); k3  = k(3); k4  = k(4);
L = l(end);

ptt = - ((i*omega).^2+(k1+k2+k3)*i*omega+k1*k3)./(vt*(i*omega+k2+k3));
pts = k1*k3./(vt*(i*omega+k2+k3));
pst = k2*k4./(vs*(i*omega+k2+k3));
pss = - ((i*omega).^2+(k2+k3+k4)*i*omega+k2*k4)./(vs*(i*omega+k2+k3));

alpha = 1/2*(ptt+pss);
beta = 1/2*sqrt((ptt-pss).^2+4*pts.*pst);

phi1 = alpha + beta;
phi2 = alpha - beta;

if vt>0 && vs>0    % parallel-flow configuration

% calculation of the frequency responses using exponential functions

Att = (phi1-pss)./(phi1-phi2);
Btt = (phi2-pss)./(phi1-phi2);
Ats = pts./(phi1-phi2);
Ast = pst./(phi1-phi2);
Ass = (phi1-ptt)./(phi1-phi2);
Bss = (phi2-ptt)./(phi1-phi2);

gttw = repmat(Att.',1,length(l)).* exp(phi1.'*l) - repmat(Btt.',1,length(l)).* exp(phi2.'*l);
gtsw = repmat(Ats.',1,length(l)).*(exp(phi1.'*l)-exp(phi2.'*l));
gstw = repmat(Ast.',1,length(l)).*(exp(phi1.'*l)-exp(phi2.'*l));
gssw = repmat(Ass.',1,length(l)).* exp(phi1.'*l) - repmat(Bss.',1,length(l)).* exp(phi2.'*l);

% calculation of the frequency responses using hyperbolic functions 

% gttw = exp(alpha.'*l).*( cosh(beta.'*l) + repmat(((alpha - pss)./beta).',1,length(l)).*sinh(beta.'*l) );
% gtsw = repmat((pts./beta).',1,length(l)) .* ( exp(alpha.'*l).*sinh(beta.'*l) );
% gstw = repmat((pst./beta).',1,length(l)) .* ( exp(alpha.'*l).*sinh(beta.'*l) );
% gssw = exp(alpha.'*l).*( cosh(beta.'*l) + repmat(((alpha - ptt)./beta).',1,length(l)).*sinh(beta.'*l) );

else 

disp('This function works only for vt>0 and vs>0.')

end
