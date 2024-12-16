function [G11h,G12h,G21h,G22h] = trans_widm_rozloz_zgodne(omega,l,LAMBDA,K)

% Funkcja obliczaj¹ca transmitancjê widmow¹ uk³adu hiperbolicznego o parametrach roz³o¿onych 
% z jedn¹ zmienn¹ przestrzenn¹ x, opisanego uk³adem dwóch liniowych równañ ró¿niczkowych postaci:
%
% dx1(l,t)/dt + lambda1*dx1(l,t)/dl = k11*x1(l,t) + k12*x2(l,t)		       
% dx2(l,t)/dt + lambda2*dx2(l,t)/dl = k21*x1(l,t) + k22*x2(l,t)		       
%
% dla warunków brzegowych okreœlonych dla l=0 oraz l=L, postaci:
% x1(0,t)=x10, x2(0,t)=x2L
%
% argumenty wejœciowe: 
% omega - wektor pulsacji, 
% l - wektor zmiennej przestrzennej,
% LAMBDA - diagonalna macierz wartoœci w³asnych (prêdkoœci) uk³adu
% K - macierz sprzêgaj¹ca
%
% wartoœci zwracane - transmitancje widmowe uk³adu:
%
% G11(l,i*omega) = X1(l,i*omega)/X1(0,i*omega)
% G12(l,i*omega) = X1(l,i*omega)/X2(L,i*omega)
% G21(l,i*omega) = X2(l,i*omega)/X1(0,i*omega)
% G22(l,i*omega) = X2(l,i*omega)/X2(L,i*omega)

lambda1 = LAMBDA(1,1);
lambda2 = LAMBDA(2,2);
k11  = K(1,1);
k12  = K(1,2);
k21  = K(2,1);
k22  = K(2,2);

%s1 = 1/2*(k11+k22+sqrt((k11-k22)^2+4*k12*k21));
%s2 = 1/2*(k11+k22-sqrt((k11-k22)^2+4*k12*k21));

%alpha = -1/(2*lambda1)*(j*omega-k11)-1/(2*lambda2)*(j*omega-k22);
%beta = 1/2*sqrt( (1/lambda1*(j*omega-k11)+ 1/lambda2*(j*omega-k22)).^2 - 4/(lambda1*lambda2)*(j*omega-s1).*(j*omega-s2) );

p11 = (k11-j*omega)/lambda1;
p12 = k12/lambda1;
p21 = k21/lambda2;
p22 = (k22-j*omega)/lambda2;

alpha = 1/2*(p11+p22);
beta = 1/2*sqrt((p11-p22).^2+4*p12*p21);

% phi1, phi22 - wartoœci w³asne macierzy P

phi1 = alpha + beta;
phi2 = alpha - beta;

A11 = (phi1-p22)./(phi1-phi2);
B11 = (phi2-p22)./(phi1-phi2);

A12 = p12./(phi1-phi2);

A21 = p21./(phi1-phi2);

A22 = (phi1-p11)./(phi1-phi2);
B22 = (phi2-p11)./(phi1-phi2);

% transmitancje operatorowe (widmowe) uk³adu  wykorzystaniem funkcji
% wyk³adniczych

G11w = repmat(A11.',1,length(l)).* exp(phi1.'*l) - repmat(B11.',1,length(l)).* exp(phi2.'*l);
G12w = repmat(A12.',1,length(l)).*(exp(phi1.'*l)-exp(phi2.'*l));
G21w = repmat(A21.',1,length(l)).*(exp(phi1.'*l)-exp(phi2.'*l));
G22w = repmat(A22.',1,length(l)).* exp(phi1.'*l) - repmat(B22.',1,length(l)).* exp(phi2.'*l);

% transmitancje operatorowe po przekszta³ceniach z wykorzystaniem funkcji hiperbolicznych

G11h = exp(alpha.'*l).*( cosh(beta.'*l) + repmat(((alpha - p22)./beta).',1,length(l)).*sinh(beta.'*l) );
G12h = repmat((p12./beta).',1,length(l)) .* ( exp(alpha.'*l).*sinh(beta.'*l) );
G21h = repmat((p21./beta).',1,length(l)) .* ( exp(alpha.'*l).*sinh(beta.'*l) );
G22h = exp(alpha.'*l).*( cosh(beta.'*l) + repmat(((alpha - p11)./beta).',1,length(l)).*sinh(beta.'*l) );
