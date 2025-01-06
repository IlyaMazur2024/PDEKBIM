
function P_s_o = P_sin_omega(G,omega,t);

% funkcja oblicza wartoœæ funkcji (P(omega)/omega)*sin(omega*t) 
% G - transmitancja widmowa (roz³o¿ona), omega - wektor pulsacji, t - czas (skalar)

[length_omega,length_x]=size(G);

P_s_o = (real(G).').*repmat(sin(omega.*t)./omega,length_x,1);
%P_s_o = ((real(G).').*sin(omega.*t))./omega



