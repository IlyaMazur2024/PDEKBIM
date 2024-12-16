function [R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_parametry(P,L,l1,Dl)

% Funkcja wyznacza warto�ci rezystancji hydraulicznej w odcinkach ruroci�gu
% o d�ugo�ci Lp z pojedynczym wyciekiem o �rednicy D1 w punkcie l1
% Ci�nienie wlotowe P

% Pozosta�e parametry ruroci�gu (d):
D = 0.2;		 % �rednica [m]
A = pi*D^2/4;	 % pole powierzchni [m^2]
alpha = 0;       % nachylenie [rad]
rho = 1000;      % g�sto�� p�ynu (woda) [kg/m^3]
ni = 10E-6;      % lepko�� kinematyczna (woda) [m^2/s]
lambda = 0.01;   % wsp�czynnik tarcia Darcy'ego - Weisbacha, 
                 % zale�ny od Re i wsp. chropowato�ci, przyj�to w. �redni�
c = 1200;        % pr�dko�� fali ci�nienia (ruroci�g stalowy) [m/s]
g = 9.81;	     % przyspieszenie ziemskie [m/s^2]

% ===========================================================================
% USZKODZENIE: w punkcie ruroci�gu l=l1 pojawia si� wyciek do atmosfery.
% W ruroci�gu nast�puje zmiana warto�ci strumieni oraz rozk�adu spadk�w ci�nie�.
l2 = L-l1;

% Wyciek do atmosfery zamodelowany jako kryza o polu powierzchni Al.
% r�wnanie kryzy: ql=kl*sqrt(dpl), gdzie: ql - strumie� masowy wycieku, 
% dpl - spadek ci�nienia na kryzie (wycieku)
% Zak�ada si� przep�yw turbulentny, Re>Recr = 10000
% przy czym Re = ql*Dl/(Al*ni*rho)

cf = 0.7;	     % �redni wsp�czynnik przep�ywu dla kryzy (wycieku)
Al = pi*Dl^2/4;  % pole powierzchni kryzy (wycieku)
kl = cf*Al*sqrt(2*rho); 

% W cz�ci odcinku przed ruroci�giem wskutek wzrostu strumienia z q do q1 
% nast�pi wzrost warto�ci spadku ci�nienia dp1
% r�wnanie spadku ci�nienia dla odcinka 1 ruroci�gu:
% q1 = k1*sqrt(dp1)

k1 = A*sqrt(2*D*rho/(lambda*l1));

% W cz�ci odcinku za ruroci�giem wskutek spadku strumienia z q do q2
% nast�pi spadek warto�ci strat ci�nienia dp2
% r�wnanie spadku ci�nienia dla odcinka 2 ruroci�gu:
% q2 = k2*sqrt(dp2)

k2 = A*sqrt(2*D*rho/(lambda*l2));

% Na zaworze wylotowym (zamodelowanym jako kryza) r�wnie� zmienia
% si� warto�� spadku ci�nienia wskutek spadku warto�ci strumienia 
% z q do q2
% R�wnanie kryzy: q2=kk*sqrt(dpk), gdzie: q2 - strumie� masowy,
% dpk - spadek ci�nienia na kryzie (zaworze)

cf = 0.7;	     % �redni wsp�czynnik przep�ywu dla kryzy (zaworu)
Dk = D/5;        % �rednica kryzy (zaworu)
Ak = pi*Dk^2/4;  % pole powierzchni kryzy (wycieku)
kk = cf*Ak*sqrt(2*rho); 

% ROZWI�ZANIE symboliczne (7 szukanych niewiadomych: dp1, dp2, dpl, dpk, q1, ql, q2)
% oraz 7 r�wna�: 4 r�wnania opisuj�ce spadek ci�nienia na elementach ruroci�gu
% 2 r�wnania bilansowe spadk�w ci�nienia oraz r�wnanie bilanosowe strumieni

syms dp1 dp2 dpl dpk q1 q2 ql positive

eq1 = 'q1 = k1*sqrt(dp1)';
eq2 = 'q2 = k2*sqrt(dp2)';
eq3 = 'ql = kl*sqrt(dpl)';
eq4 = 'q2 = kk*sqrt(dpk)';
eq5 = 'dp1+dpl = P';
eq6 = 'dp1+dp2+dpk = P';
eq7 = 'q1 = q2+ql';

[dp1,dp2,dpl,dpk,q1,q2,ql] = solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,'dp1,dp2,dpl,dpk,q1,q2,ql');
dp1 = subs(dp1(4));
dp2 = subs(dp2(4));
dpl = subs(dpl(4));
dpk = subs(dpk(4));
q1 = subs(q1(4));
q2 = subs(q2(4));
ql = subs(ql(4));


Rl = ql/(2*rho*cf^2*Al^2)       % rezystancja hydrauliczna wycieku 
Rl = dpl/ql
Rel = ql*Dl/(Al*ni*rho)  	     % liczba Reynoldsa wycieku

R1 = lambda*q1*l1/(2*rho*D*A^2) % rezystancja odcinka 1 ruroci�gu
R1 = dp1/q1
Re1 = q1*D/(A*ni*rho)	         % liczba Reynoldsa odcinka 1 ruroci�gu

R2 = lambda*q2*l2/(2*rho*D*A^2) % rezystancja odcinka 2 ruroci�gu
R2 = dp2/q2
Re2 = q2*D/(A*ni*rho)	         % liczba Reynoldsa odcinka 2 ruroci�gu

Rk = qk/(2*rho*cf^2*Ak^2)       % rezystancja hydrauliczna zaworu ko�cowego 
Rk = dpk/q2
Rek = q2*Dk/(Ak*ni*rho)	      % liczba Reynoldsa zaworu ko�cowego
