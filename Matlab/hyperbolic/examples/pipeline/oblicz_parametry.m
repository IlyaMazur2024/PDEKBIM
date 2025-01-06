function [R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_parametry(P,L,l1,Dl)

% Funkcja wyznacza wartoœci rezystancji hydraulicznej w odcinkach ruroci¹gu
% o d³ugoœci Lp z pojedynczym wyciekiem o œrednicy D1 w punkcie l1
% Ciœnienie wlotowe P

% Pozosta³e parametry ruroci¹gu (d):
D = 0.2;		 % œrednica [m]
A = pi*D^2/4;	 % pole powierzchni [m^2]
alpha = 0;       % nachylenie [rad]
rho = 1000;      % gêstoœæ p³ynu (woda) [kg/m^3]
ni = 10E-6;      % lepkoœæ kinematyczna (woda) [m^2/s]
lambda = 0.01;   % wspó³czynnik tarcia Darcy'ego - Weisbacha, 
                 % zale¿ny od Re i wsp. chropowatoœci, przyjêto w. œredni¹
c = 1200;        % prêdkoœæ fali ciœnienia (ruroci¹g stalowy) [m/s]
g = 9.81;	     % przyspieszenie ziemskie [m/s^2]

% ===========================================================================
% USZKODZENIE: w punkcie ruroci¹gu l=l1 pojawia siê wyciek do atmosfery.
% W ruroci¹gu nastêpuje zmiana wartoœci strumieni oraz rozk³adu spadków ciœnieñ.
l2 = L-l1;

% Wyciek do atmosfery zamodelowany jako kryza o polu powierzchni Al.
% równanie kryzy: ql=kl*sqrt(dpl), gdzie: ql - strumieñ masowy wycieku, 
% dpl - spadek ciœnienia na kryzie (wycieku)
% Zak³ada siê przep³yw turbulentny, Re>Recr = 10000
% przy czym Re = ql*Dl/(Al*ni*rho)

cf = 0.7;	     % œredni wspó³czynnik przep³ywu dla kryzy (wycieku)
Al = pi*Dl^2/4;  % pole powierzchni kryzy (wycieku)
kl = cf*Al*sqrt(2*rho); 

% W czêœci odcinku przed ruroci¹giem wskutek wzrostu strumienia z q do q1 
% nast¹pi wzrost wartoœci spadku ciœnienia dp1
% równanie spadku ciœnienia dla odcinka 1 ruroci¹gu:
% q1 = k1*sqrt(dp1)

k1 = A*sqrt(2*D*rho/(lambda*l1));

% W czêœci odcinku za ruroci¹giem wskutek spadku strumienia z q do q2
% nast¹pi spadek wartoœci strat ciœnienia dp2
% równanie spadku ciœnienia dla odcinka 2 ruroci¹gu:
% q2 = k2*sqrt(dp2)

k2 = A*sqrt(2*D*rho/(lambda*l2));

% Na zaworze wylotowym (zamodelowanym jako kryza) równie¿ zmienia
% siê wartoœæ spadku ciœnienia wskutek spadku wartoœci strumienia 
% z q do q2
% Równanie kryzy: q2=kk*sqrt(dpk), gdzie: q2 - strumieñ masowy,
% dpk - spadek ciœnienia na kryzie (zaworze)

cf = 0.7;	     % œredni wspó³czynnik przep³ywu dla kryzy (zaworu)
Dk = D/5;        % œrednica kryzy (zaworu)
Ak = pi*Dk^2/4;  % pole powierzchni kryzy (wycieku)
kk = cf*Ak*sqrt(2*rho); 

% ROZWI¥ZANIE symboliczne (7 szukanych niewiadomych: dp1, dp2, dpl, dpk, q1, ql, q2)
% oraz 7 równañ: 4 równania opisuj¹ce spadek ciœnienia na elementach ruroci¹gu
% 2 równania bilansowe spadków ciœnienia oraz równanie bilanosowe strumieni

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

R1 = lambda*q1*l1/(2*rho*D*A^2) % rezystancja odcinka 1 ruroci¹gu
R1 = dp1/q1
Re1 = q1*D/(A*ni*rho)	         % liczba Reynoldsa odcinka 1 ruroci¹gu

R2 = lambda*q2*l2/(2*rho*D*A^2) % rezystancja odcinka 2 ruroci¹gu
R2 = dp2/q2
Re2 = q2*D/(A*ni*rho)	         % liczba Reynoldsa odcinka 2 ruroci¹gu

Rk = qk/(2*rho*cf^2*Ak^2)       % rezystancja hydrauliczna zaworu koñcowego 
Rk = dpk/q2
Rek = q2*Dk/(Ak*ni*rho)	      % liczba Reynoldsa zaworu koñcowego
