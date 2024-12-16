
function [G] = trans_widm(omega,l,q)

% Funkcja obliczaj¹ca transmitancjê widmow¹ modelu ruroci¹gu w punkcie l
% dla danej wartoœci strumienia masowego q
% G = -Zc*( (exp(gamma*(Lp+l)) - exp(gamma*(Lp-l))) / (exp(2*gamma*Lp) + 1) )

% Parametry ruroci¹gu:
Lp=1000;     % d³ugoœæ ruroci¹gu
D=0.2;		 % œrednica ruroci¹gu
A=pi*D^2/4;  % pole powierzchni ruroci¹gu
rho=1000;    % gêstoœæ p³ynu (stal)
lambda=0.01; % wspó³czynnik tarcia Darcy'ego - Weisbacha
c=1200;      % prêdkoœæ fali ciœnienia (ruroci¹g stalowy)
g=9.81;	    % przyspieszenie ziemskie

R=lambda*q/(rho*D*A^2); 
L=1/A; 					
C=A/c^2;

Zc=sqrt((R+j*omega*L)./(j*omega*C));			% impedancja charakterystyczna ruroci¹gu
gamma = sqrt((R+j*omega*L).*(j*omega*C)); % sta³a rozchodzenia siê fal


%G = -repmat(Zc.',1,length(l)).*tanh(gamma.'*l);

Ga = -repmat(Zc.',1,length(l)).*(exp(gamma.'*(Lp+l)) - exp(gamma.'*(Lp-l)));
Gb = repmat(exp(2*gamma.'*Lp)+1,1,length(l));

G = Ga ./ Gb;