
function [G] = trans_widm(omega,l,q)

% Funkcja obliczaj�ca transmitancj� widmow� modelu ruroci�gu w punkcie l
% dla danej warto�ci strumienia masowego q
% G = -Zc*( (exp(gamma*(Lp+l)) - exp(gamma*(Lp-l))) / (exp(2*gamma*Lp) + 1) )

% Parametry ruroci�gu:
Lp=1000;     % d�ugo�� ruroci�gu
D=0.2;		 % �rednica ruroci�gu
A=pi*D^2/4;  % pole powierzchni ruroci�gu
rho=1000;    % g�sto�� p�ynu (stal)
lambda=0.01; % wsp�czynnik tarcia Darcy'ego - Weisbacha
c=1200;      % pr�dko�� fali ci�nienia (ruroci�g stalowy)
g=9.81;	    % przyspieszenie ziemskie

R=lambda*q/(rho*D*A^2); 
L=1/A; 					
C=A/c^2;

Zc=sqrt((R+j*omega*L)./(j*omega*C));			% impedancja charakterystyczna ruroci�gu
gamma = sqrt((R+j*omega*L).*(j*omega*C)); % sta�a rozchodzenia si� fal


%G = -repmat(Zc.',1,length(l)).*tanh(gamma.'*l);

Ga = -repmat(Zc.',1,length(l)).*(exp(gamma.'*(Lp+l)) - exp(gamma.'*(Lp-l)));
Gb = repmat(exp(2*gamma.'*Lp)+1,1,length(l));

G = Ga ./ Gb;