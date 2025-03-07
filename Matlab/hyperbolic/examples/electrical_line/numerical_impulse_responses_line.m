
function [g11n,g12n,g21n,g22n] = numerical_impulse_responses_line(t,l,R,L,G,C)
    
% =========================================================================
% 
% Function m-file name: numerical_impulse_responses.m
%
% Purpose: function calculates spatiotemporal impulse reponses of the 2x2 
% weakly coupled hyperbolic system based on the numerical calculation
% of the Bromwich integral representing the inverse Laplace transform.
% 
% Syntax: [g11,g12,g21,g22]=numerical_impulse_responses(LAMBDA,K,t,l)
%
% Function parameters:
% LAMBDA - 2 x 2 matrix of eigenvalues,
% K - 2 x 2 coupling matrix, 
% t - vector of time instants,
% l - vector of spatial positions,
% L - length of the spatial domain.
%
% =========================================================================

g11n = zeros(length(t),length(l));  
g12n = zeros(length(t),length(l));  
g21n = zeros(length(t),length(l));  
g22n = zeros(length(t),length(l));  
                                     
% the integration line Re(s)=x, (x-iy,x+iy)
x = 0;
y = 1000;
step = 0.001;                  % step along the integration line 
tau = [0:step:1];

s_tau = x + 1i*(-y+2*y*tau);  % parametrization of the integration line
s_tau_deriv = 2*y*1i;               

[dum,dum,dum,dum,G11,G12,G21,G22] = freq_resp_trans_line(s_tau/1i,l,R,L,G,C);

for j=1:length(t)                   %  successive time steps
  
% transfer function G_11(l,s)
  G11_tau = G11.*repmat(exp(s_tau.'*t(j)),1,length(l))*s_tau_deriv;
  int_G11 = sum(G11_tau)*step;
  g11n(j,:) = 1/(2*pi*1i)*int_G11;
  
  % transfer function G_12(l,s) 
  G12_tau = G12.*repmat(exp(s_tau.'*t(j)),1,length(l))*s_tau_deriv;
  int_G12 = sum(G12_tau)*step;
  g12n(j,:) = 1/(2*pi*1i)*int_G12;
  
% transfer function G_21(l,s) 
  G21_tau = G21.*repmat(exp(s_tau.'*t(j)),1,length(l))*s_tau_deriv;
  int_G21 = sum(G21_tau)*step;
  g21n(j,:) = 1/(2*pi*1i)*int_G21;  

% transfer function G_22(l,s)
  G22_tau = G22.*repmat(exp(s_tau.'*t(j)),1,length(l))*s_tau_deriv;
  int_G22 = sum(G22_tau)*step;
  g22n(j,:) = 1/(2*pi*1i)*int_G22;
  
end
