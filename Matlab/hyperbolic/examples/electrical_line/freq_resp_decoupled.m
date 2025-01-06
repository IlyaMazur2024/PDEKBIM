function [G11c,G12c,G21c,G22c,G11i,G12i,G21i,G22i] = freq_resp_decoupled(omega,l,LAMBDA,K)

% Funkcja obliczaj¹ca transmitancjê widmow¹ uk³adu hiperbolicznego o parametrach roz³o¿onych 
% z jedn¹ zmienn¹ przestrzenn¹ x, opisanego uk³adem dwóch liniowych równañ ró¿niczkowych postaci:
% 
% dx1(l,t)/dt + lambda1*dx1(l,t)/dl = k11*x1(l,t) + k12*x2(l,t)		       
% dx2(l,t)/dt + lambda2*dx2(l,t)/dl = k21*x1(l,t) + k22*x2(l,t)		       
%
% dla warunków brzegowych okreœlonych dla x=0, postaci: x1(0,t)=x10, x2(0,t)=x20
%
% argumenty wejœciowe: 
% omega - wektor pulsacji, 
% l - wektor zmiennej przestrzennej,
%
% wartoœci zwracane - transmitancje operatorowe uk³adu:
%
% G11(l,s) = X1(l,s)/X1(0,s)
% G12(l,s) = X1(l,s)/X2(0,s)
% G21(l,s) = X2(l,s)/X1(0,s)
% G22(l,s) = X2(l,s)/X2(0,s)

L = l(end);

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


% =================================================================================
% warunki brzegowe zgodne

A11 = (phi1-p22)./(phi1-phi2);
B11 = (phi2-p22)./(phi1-phi2);
A12 = p12./(phi1-phi2);
A21 = p21./(phi1-phi2);
A22 = (phi1-p11)./(phi1-phi2);
B22 = (phi2-p11)./(phi1-phi2);

% transmitancje operatorowe (widmowe) uk³adu z wykorzystaniem funkcji wyk³adniczych

% poni¿sze dzia³aj¹ tylko dla skalarnego l
% G11c = A11 .* exp(q1*l) - B11 .* exp(q2*l);  
% G12c = A12.*(exp(q1.*l)-exp(q2.*l));       
% G21c = A21.*(exp(q1.*l)-exp(q2.*l));       
% G22c = A2.*exp(q1*l) - B2.* exp(q2*l);    

% poni¿sze dzia³aj¹ równie¿ dla wektorowego l
G11c = repmat(A11.',1,length(l)).* exp(phi1.'*l) - repmat(B11.',1,length(l)).* exp(phi2.'*l);
G12c = repmat(A12.',1,length(l)).*(exp(phi1.'*l)-exp(phi2.'*l));
G21c = repmat(A21.',1,length(l)).*(exp(phi1.'*l)-exp(phi2.'*l));
G22c = repmat(A22.',1,length(l)).* exp(phi1.'*l) - repmat(B22.',1,length(l)).* exp(phi2.'*l);

% transmitancje operatorowe po przekszta³ceniach z wykorzystaniem funkcji hiperbolicznych

% poni¿sze dzia³aj¹ tylko dla skalarnego l
% G11hi = exp(alpha*l).*( cosh(beta*l) + ((alpha - p22)./beta).*sinh(beta*l)); 
% G12hi = p12./beta.*exp(alpha*l).*sinh(beta*l);  
% G21h = p21./beta.*exp(alpha.*l).*sinh(beta.*l); 
% G22h = exp(alpha.*l).*( cosh(beta.*l) + ((alpha - p11)./beta).*sinh(beta.*l) ); 

% poni¿sze dzia³aj¹ równie¿ dla wektorowego l
% G11hi = exp(alpha.'*l).*( cosh(beta.'*l) + repmat(((alpha - p22)./beta).',1,length(l)).*sinh(beta.'*l) );
% G12hi = repmat((p12./beta).',1,length(l)) .* ( exp(alpha.'*l).*sinh(beta.'*l) );
% G21hi = repmat((p21./beta).',1,length(l)) .* ( exp(alpha.'*l).*sinh(beta.'*l) );
% G22hi = exp(alpha.'*l).*( cosh(beta.'*l) + repmat(((alpha - p11)./beta).',1,length(l)).*sinh(beta.'*l) );


% =================================================================================
% warunki brzegowe przeciwne

A11 = (exp(phi2*L).*(phi1-p22)) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));
B11 = (exp(phi1*L).*(phi2-p22)) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));
A12 = p12 ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));
A21 = p21*exp(phi2*L) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));
B21 = p21*exp(phi1*L) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));
A22 = (phi2-p11) ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));
B22 = (phi1-p11) ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));

% transmitancje operatorowe (widmowe) uk³adu  wykorzystaniem funkcji
% wyk³adniczych

% poni¿sze dzia³aj¹ tylko dla skalarnego l
% G11i = A11 .* exp(q1*l) - B11 .* exp(q2*l);  
% G12i = A12.*(exp(q2.*l)-exp(q1.*l));       
% G21i = A21.* exp(q1.*l) - B21 .* exp(q2.*l) ;       
% G22i = A2.*exp(q2*l) - B2.* exp(q1*l);    

% poni¿sze dzia³aj¹ równie¿ dla wektorowego l
G11i = repmat(A11.',1,length(l)).* exp(phi1.'*l) - repmat(B11.',1,length(l)).* exp(phi2.'*l);
G12i = repmat(A12.',1,length((l))).*(exp(phi2.'*(l))-exp(phi1.'*(l)));
G21i = repmat(A21.',1,length(l)).* exp(phi1.'*l) - repmat(B21.',1,length(l)).* exp(phi2.'*l);
G22i = (repmat(A22.',1,length(l)).* exp(phi2.'*(l)) - repmat(B22.',1,length((l))).* exp(phi1.'*(l)));


% transmitancje operatorowe po przekszta³ceniach z wykorzystaniem funkcji hiperbolicznych

C11 = beta./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );
D11 = (alpha-p22)./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );
C12 = p12./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );
C21 = p21./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );
C22 = beta./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );
D22 = (alpha-p11)./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );

% poni¿sze dzia³aj¹ tylko dla skalarnego l
% G11hi =   ( exp(alpha*l).* ( beta.*cosh(beta*(l-L)) + (alpha-p22).*sinh(beta*(l-L))) ) ./ ...
%         ( beta .* cosh(beta*L)-(alpha-p22).*sinh(beta*L) ); 
% G12hi =   ( p12*exp(alpha*(l-L)).*sinh(beta*l) ) ./ ... 
%         ( beta.*cosh(beta*L) + (alpha-p11).*sinh(beta*L) ); 
% G21hi =   ( p21*exp(alpha.*l).*sinh(beta.*(l-L)) ) ./ ... 
%         ( beta.*cosh(beta*L) - (alpha-p22).*sinh(beta*L) ); 
% G22hi =   ( exp(alpha*(l-L)).* ( beta.*cosh(beta*l) + (alpha-p11).*sinh(beta*l)) ) ./ ...
%         ( beta .* cosh(beta*L)+(alpha-p11).*sinh(beta*L) ); dzia³a tylko dla skalarnego l

% poni¿sze dzia³aj¹ równie¿ dla wektorowego l
% G11hi = repmat(C11.',1,length(l)) .* (exp(alpha.'*(l)).*cosh(beta.'*(l-L))) + ...
%       repmat(D11.',1,length(l)) .* (exp(alpha.'*(l)).*sinh(beta.'*(l-L))); 
% G12hi = repmat(C12.',1,length(l)) .* (exp(alpha.'*(l-L)).*sinh(beta.'*l));
% G21hi = repmat(C21.',1,length(l)) .* (exp(alpha.'*(l)).*sinh(beta.'*(l-L)));
% G22hi = repmat(C22.',1,length(l)) .* (exp(alpha.'*(l-L)).*cosh(beta.'*l)) + ...
%       repmat(D22.',1,length(l)) .* (exp(alpha.'*(l-L)).*sinh(beta.'*l)); 

