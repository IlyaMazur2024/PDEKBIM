function [G11w,G12w,G21w,G22w] = calc_hyper_distr_freq_tf(omega,l,L,LAMBDA,K)

% Funkcja obliczaj¹ca transmitancjê widmow¹ uk³adu hiperbolicznego o parametrach roz³o¿onych 
% z jedn¹ zmienn¹ przestrzenn¹ x, opisanego uk³adem dwóch liniowych równañ ró¿niczkowych postaci:
%
% dx1(l,t)/dt + lambda1*dx1(l,t)/dl = k11*x1(l,t) + k12*x2(l,t)		       
% dx2(l,t)/dt + lambda2*dx2(l,t)/dl = k21*x1(l,t) + k22*x2(l,t)		       
%
% dla warunków brzegowych okreœlonych dla lambda1>0 i lambda2>0 
% w postaci: x1(0,t)=x10, x2(0,t)=x20 oraz 
% dla warunków brzegowych okreœlonych dla lambda1>0 i lambda2<0 
% w postaci: x1(0,t)=x10, x2(L,t)=x2L
%
% argumenty wejœciowe: 
% omega - wektor zawieraj¹cy wartoœci pulsacji, 
% l - wektor zawieraj¹cy wartoœci zmiennej przestrzennej,
% L - d³ugoœæ uk³adu,
% LAMBDA - diagonalna macierz wartoœci w³asnych (prêdkoœci) uk³adu,
% K - macierz sprzêgaj¹ca.
%
% wartoœci zwracane - transmitancje widmowe uk³adu:
%
% G11(l,i*omega) = X1(l,i*omega)/X1(0,i*omega)
% G12(l,i*omega) = X1(l,i*omega)/X2(0,i*omega)
% G21(l,i*omega) = X2(l,i*omega)/X1(0,i*omega)
% G22(l,i*omega) = X2(l,i*omega)/X2(0,i*omega)

lambda1 = LAMBDA(1,1);
lambda2 = LAMBDA(2,2);
k11  = K(1,1);
k12  = K(1,2);
k21  = K(2,1);
k22  = K(2,2);

p11 = (k11-j*omega)/lambda1;
p12 = k12/lambda1;
p21 = k21/lambda2;
p22 = (k22-j*omega)/lambda2;

alpha = 1/2*(p11+p22);
beta = 1/2*sqrt((p11-p22).^2+4*p12*p21);

% phi1, phi22 - wartoœci w³asne macierzy P

phi1 = alpha + beta;
phi2 = alpha - beta;

if lambda1>0 && lambda2>0 % collocated boundary conditions  

A11 = (phi1-p22)./(phi1-phi2);
B11 = (phi2-p22)./(phi1-phi2);

A12 = p12./(phi1-phi2);

A21 = p21./(phi1-phi2);

A22 = (phi1-p11)./(phi1-phi2);
B22 = (phi2-p11)./(phi1-phi2);

% transmitancje operatorowe (widmowe) uk³adu z wykorzystaniem funkcji wyk³adniczych

G11w = repmat(A11.',1,length(l)).* exp(phi1.'*l) - repmat(B11.',1,length(l)).* exp(phi2.'*l);
G12w = repmat(A12.',1,length(l)).*(exp(phi1.'*l)-exp(phi2.'*l));
G21w = repmat(A21.',1,length(l)).*(exp(phi1.'*l)-exp(phi2.'*l));
G22w = repmat(A22.',1,length(l)).* exp(phi1.'*l) - repmat(B22.',1,length(l)).* exp(phi2.'*l);

% transmitancje operatorowe po przekszta³ceniach z wykorzystaniem funkcji hiperbolicznych

%G11w = exp(alpha.'*l).*( cosh(beta.'*l) + repmat(((alpha - p22)./beta).',1,length(l)).*sinh(beta.'*l) );
%G12w = repmat((p12./beta).',1,length(l)) .* ( exp(alpha.'*l).*sinh(beta.'*l) );
%G21w = repmat((p21./beta).',1,length(l)) .* ( exp(alpha.'*l).*sinh(beta.'*l) );
%G22w = exp(alpha.'*l).*( cosh(beta.'*l) + repmat(((alpha - p11)./beta).',1,length(l)).*sinh(beta.'*l) );

elseif lambda1>0 && lambda2<0 % anti-collocated boundary conditions
    
A11 = (exp(phi2*L).*(phi1-p22)) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));
B11 = (exp(phi1*L).*(phi2-p22)) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));

A12 = p12 ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));

A21 = p21*exp(phi2*L) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));
B21 = p21*exp(phi1*L) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));

A22 = (phi2-p11) ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));
B22 = (phi1-p11) ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));

% transmitancje operatorowe (widmowe) uk³adu  wykorzystaniem funkcji
% wyk³adniczych

G11w = repmat(A11.',1,length(l)).* exp(phi1.'*l) - repmat(B11.',1,length(l)).* exp(phi2.'*l);
G12w = repmat(A12.',1,length((l))).*(exp(phi2.'*(l))-exp(phi1.'*(l)));
G21w = repmat(A21.',1,length(l)).* exp(phi1.'*l) - repmat(B21.',1,length(l)).* exp(phi2.'*l);
G22w = (repmat(A22.',1,length(l)).* exp(phi2.'*(l)) - repmat(B22.',1,length((l))).* exp(phi1.'*(l)));

% transmitancje operatorowe po przekszta³ceniach z wykorzystaniem funkcji hiperbolicznych

% C11 = beta./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );
% D11 = (alpha-p22)./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );
% C12 = p12./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );
% C21 = p21./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );
% C22 = beta./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );
% D22 = (alpha-p11)./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );
% 
% G11w = repmat(C11.',1,length(l)) .* (exp(alpha.'*(l)).*cosh(beta.'*(l-L))) + ...
%        repmat(D11.',1,length(l)) .* (exp(alpha.'*(l)).*sinh(beta.'*(l-L))); 
% G12w = repmat(C12.',1,length(l)) .* (exp(alpha.'*(l-L)).*sinh(beta.'*l));
% G21w = repmat(C21.',1,length(l)) .* (exp(alpha.'*(l)).*sinh(beta.'*(l-L)));
% G22w = repmat(C22.',1,length(l)) .* (exp(alpha.'*(l-L)).*cosh(beta.'*l)) + ...
%        repmat(D22.',1,length(l)) .* (exp(alpha.'*(l-L)).*sinh(beta.'*l)); 

end
