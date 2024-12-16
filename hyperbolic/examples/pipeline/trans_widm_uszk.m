
function [Gvl] = trans_widm_uszk(omega,l1,Dl,q1,q2,ql)

% Funkcja obliczaj¹ca transmitancjê widmow¹ modelu ruroci¹gu (dla punktu koñcowego l=Lp)
% l1 - po³o¿enie uszkodzenia od pocz¹tku ruroci¹gu, 0<l1<Lp
% Dl - œrednica wycieku  
% q1 - strumieñ masowy cieczy w odcinku przed wyciekiem
% q2 - strumieñ masowy cieczy w odcinku za wyciekiem
% ql - strumieñ masowy wycieku

% Parametry ruroci¹gu:
%
Lp=1000;     % d³ugoœæ ruroci¹gu
l2=Lp-l1;
D=0.2;		 % œrednica ruroci¹gu
A=pi*D^2/4;  % pole powierzchni ruroci¹gu
rho=1000;    % gêstoœæ p³ynu (stal)
lambda=0.01; % wspó³czynnik tarcia Darcy'ego - Weisbacha
c=1200;      % prêdkoœæ fali ciœnienia (ruroci¹g stalowy)
g=9.81;	    % przyspieszenie ziemskie
q=27.29;	    % strumieñ masowy cieczy (ustalony punkt pracy)

R1=lambda*q1/(rho*D*A^2); 
R2=lambda*q2/(rho*D*A^2); 
L=1/A; 					
C=A/c^2;

Zc1=sqrt((R1+j*omega*L)./(j*omega*C));			% impedancja charakterystyczna odcinka ruroci¹gu przed wyciekiem
gamma1 = sqrt((R1+j*omega*L).*(j*omega*C));  % sta³a rozchodzenia siê fal w odcinku ruroci¹gu przed wyciekiem

Zc2=sqrt((R2+j*omega*L)./(j*omega*C));			% impedancja charakterystyczna odcinka ruroci¹gu za wyciekiem
gamma2 = sqrt((R2+j*omega*L).*(j*omega*C));  % sta³a rozchodzenia siê fal w odcinku ruroci¹gu za wyciekiem

cf = 0.7;	                % œredni wspó³czynnik przep³ywu dla kryzy (wycieku)
Al = pi*Dl^2/4;			    % pole powierzchni wycieku
Rl = ql/(2*rho*cf^2*Al^2);  % rezystancja hydrauliczna wycieku 
Zl = Rl;							 % wyciek ma charakter "rezystancyjny"


Gvl = ( -Zl.*Zc1.*tanh(gamma1*l1)-Zl.*Zc2.*tanh(gamma2*l2)-Zc1.*Zc2.*tanh(gamma1*l1).*tanh(gamma2*l2) ) ...
       ./ (Zl+Zc1.*tanh(gamma1*l1)+Zl.*(Zc1./Zc2).*tanh(gamma1*l1).*tanh(gamma2*l2));
Gvl = Gvl';

