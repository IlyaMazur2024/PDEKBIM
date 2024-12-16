
function [Gvl] = trans_widm_uszk(omega,l1,Dl,q1,q2,ql)

% Funkcja obliczaj�ca transmitancj� widmow� modelu ruroci�gu (dla punktu ko�cowego l=Lp)
% l1 - po�o�enie uszkodzenia od pocz�tku ruroci�gu, 0<l1<Lp
% Dl - �rednica wycieku  
% q1 - strumie� masowy cieczy w odcinku przed wyciekiem
% q2 - strumie� masowy cieczy w odcinku za wyciekiem
% ql - strumie� masowy wycieku

% Parametry ruroci�gu:
%
Lp=1000;     % d�ugo�� ruroci�gu
l2=Lp-l1;
D=0.2;		 % �rednica ruroci�gu
A=pi*D^2/4;  % pole powierzchni ruroci�gu
rho=1000;    % g�sto�� p�ynu (stal)
lambda=0.01; % wsp�czynnik tarcia Darcy'ego - Weisbacha
c=1200;      % pr�dko�� fali ci�nienia (ruroci�g stalowy)
g=9.81;	    % przyspieszenie ziemskie
q=27.29;	    % strumie� masowy cieczy (ustalony punkt pracy)

R1=lambda*q1/(rho*D*A^2); 
R2=lambda*q2/(rho*D*A^2); 
L=1/A; 					
C=A/c^2;

Zc1=sqrt((R1+j*omega*L)./(j*omega*C));			% impedancja charakterystyczna odcinka ruroci�gu przed wyciekiem
gamma1 = sqrt((R1+j*omega*L).*(j*omega*C));  % sta�a rozchodzenia si� fal w odcinku ruroci�gu przed wyciekiem

Zc2=sqrt((R2+j*omega*L)./(j*omega*C));			% impedancja charakterystyczna odcinka ruroci�gu za wyciekiem
gamma2 = sqrt((R2+j*omega*L).*(j*omega*C));  % sta�a rozchodzenia si� fal w odcinku ruroci�gu za wyciekiem

cf = 0.7;	                % �redni wsp�czynnik przep�ywu dla kryzy (wycieku)
Al = pi*Dl^2/4;			    % pole powierzchni wycieku
Rl = ql/(2*rho*cf^2*Al^2);  % rezystancja hydrauliczna wycieku 
Zl = Rl;							 % wyciek ma charakter "rezystancyjny"


Gvl = ( -Zl.*Zc1.*tanh(gamma1*l1)-Zl.*Zc2.*tanh(gamma2*l2)-Zc1.*Zc2.*tanh(gamma1*l1).*tanh(gamma2*l2) ) ...
       ./ (Zl+Zc1.*tanh(gamma1*l1)+Zl.*(Zc1./Zc2).*tanh(gamma1*l1).*tanh(gamma2*l2));
Gvl = Gvl';

