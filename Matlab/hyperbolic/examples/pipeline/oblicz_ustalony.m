function [ddp1,ddp2,ddpl,ddpk,dq1,dq2,dql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,L,l1,Dl)

% Funkcja wyznacza rozk³ad spadków ciœnieñ oraz przep³ywów wzd³u¿ osi ruroci¹gu
% o d³ugoœci Lp z pojedynczym wyciekiem o œrednicy D1 w punkcie l1
% Ciœnienie wlotowe P

% Pozosta³e parametry ruroci¹gu (d):
D=0.2;		 % œrednica [m]
A=pi*D^2/4;	 % pole powierzchni [m^2]
alpha=0;     % nachylenie [rad]
rho=1000;    % gêstoœæ p³ynu (woda) [kg/m^3]
ni=10-6;     % lepkoœæ kinematyczna (woda) [m^2/s]
lambda=0.01; % wspó³czynnik tarcia Darcy'ego - Weisbacha, 
             % zale¿ny od Re i wsp. chropowatoœci, przyjêto w. œredni¹
c=1200;      % prêdkoœæ fali ciœnienia (ruroci¹g stalowy) [m/s]
g=9.81;	    % przyspieszenie ziemskie [m/s^2]

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
%Dk = D/2;
Ak = pi*Dk^2/4;  % pole powierzchni kryzy 
kk = cf*Ak*sqrt(2*rho); 

% ROZWI¥ZANIE symboliczne (7 szukanych niewiadomych: dp1, dp2, dpl, dpk, q1, ql, q2)
% oraz 7 równañ: 4 równania opisuj¹ce spadek ciœnienia na elementach ruroci¹gu
% 2 równania bilansowe spadków ciœnienia oraz równanie bilanosowe strumieni

syms dp1 dp2 dpl dpk q1 q2 ql 

eq1 = 'q1 = k1*sqrt(dp1)';
eq2 = 'q2 = k2*sqrt(dp2)';
eq3 = 'ql = kl*sqrt(dpl)';
eq4 = 'q2 = kk*sqrt(dpk)';
eq5 = 'dp1+dpl = P';
eq6 = 'dp1+dp2+dpk = P';
eq7 = 'q1 = q2+ql';

[dp1,dp2,dpl,dpk,q1,q2,ql] = solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,'dp1,dp2,dpl,dpk,q1,q2,ql');
ddp1 = subs(dp1(1));
ddp2 = subs(dp2(1));
ddpl = subs(dpl(1));
ddpk = subs(dpk(1));
dq1 = subs(q1(1));
dq2 = subs(q2(1));
dql = subs(ql(1));

Rl = dql/(2*rho*cf^2*Al^2);       % rezystancja hydrauliczna wycieku 
Rel = dql*Dl/(Al*ni*rho);  	      % liczba Reynoldsa wycieku

R1 = lambda*dq1*l1/(2*rho*D*A^2); % rezystancja odcinka 1 ruroci¹gu
Re1 = dq1*D/(A*ni*rho);	         % liczba Reynoldsa odcinka 1 ruroci¹gu

R2 = lambda*dq2*l2/(2*rho*D*A^2); % rezystancja odcinka 2 ruroci¹gu
Re2 = dq2*D/(A*ni*rho);	         % liczba Reynoldsa odcinka 2 ruroci¹gu

Rk = dq2/(2*rho*cf^2*Ak^2);       % rezystancja hydrauliczna zaworu koñcowego 
Rek = dq2*Dk/(Ak*ni*rho);	      % liczba Reynoldsa zaworu koñcowego

