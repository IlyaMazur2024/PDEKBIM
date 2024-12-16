function [G11h,G12h,G21h,G22h] = trans_widm_rozloz_przeciwne(omega,l,L,LAMBDA,K)

% Funkcja obliczaj�ca transmitancj� widmow� uk�adu hiperbolicznego o parametrach roz�o�onych 
% z jedn� zmienn� przestrzenn� x, opisanego uk�adem dw�ch liniowych r�wna� r�niczkowych postaci:
%
% dx1(l,t)/dt + lambda1*dx1(l,t)/dl = k11*x1(l,t) + k12*x2(l,t)		       
% dx2(l,t)/dt + lambda2*dx2(l,t)/dl = k21*x1(l,t) + k22*x2(l,t)		       
%
% dla warunk�w brzegowych okre�lonych dla l=0 oraz l=L, postaci:
% x1(0,t)=x10, x2(L,t)=x2L
%
% argumenty wej�ciowe: 
% omega - wektor pulsacji, 
% l - wektor zmiennej przestrzennej,
% LAMBDA - diagonalna macierz warto�ci w�asnych (pr�dko�ci) uk�adu
% K - macierz sprz�gaj�ca
%
% warto�ci zwracane - transmitancje widmowe uk�adu:
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

% phi1, phi22 - warto�ci w�asne macierzy P

phi1 = alpha + beta;
phi2 = alpha - beta;

A11 = (exp(phi2*L).*(phi1-p22)) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));
B11 = (exp(phi1*L).*(phi2-p22)) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));

A12 = p12 ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));

A21 = p21*exp(phi2*L) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));
B21 = p21*exp(phi1*L) ./ (exp(phi2*L).*(phi1-p22) - exp(phi1*L).*(phi2-p22));

A22 = (phi2-p11) ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));
B22 = (phi1-p11) ./ (exp(phi2*L).*(phi2-p11) - exp(phi1*L).*(phi1-p11));

% transmitancje operatorowe (widmowe) uk�adu  wykorzystaniem funkcji
% wyk�adniczych

% G11 = A11 .* exp(q1*l) - B11 .* exp(q2*l);  dzia�a tylko dla skalarnego l
G11 = repmat(A11.',1,length(l)).* exp(phi1.'*l) - repmat(B11.',1,length(l)).* exp(phi2.'*l);

% G12 = A12.*(exp(q2.*l)-exp(q1.*l));       dzia�a tylko dla skalarnego l
G12 = repmat(A12.',1,length((l))).*(exp(phi2.'*(l))-exp(phi1.'*(l)));

% G21 = A21.* exp(q1.*l) - B21 .* exp(q2.*l) ;       dzia�a tylko dla skalarnego l
G21 = repmat(A21.',1,length(l)).* exp(phi1.'*l) - repmat(B21.',1,length(l)).* exp(phi2.'*l);

% G22 = A2.*exp(q2*l) - B2.* exp(q1*l);    dzia�a tylko dla skalarnego l
G22 = (repmat(A22.',1,length(l)).* exp(phi2.'*(l)) - repmat(B22.',1,length((l))).* exp(phi1.'*(l)));


% transmitancje operatorowe po przekszta�ceniach z wykorzystaniem funkcji hiperbolicznych

C11 = beta./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );
D11 = (alpha-p22)./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );

C12 = p12./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );

C21 = p21./( beta.*cosh(beta*L)-(alpha-p22).*sinh(beta*L) );

C22 = beta./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );
D22 = (alpha-p11)./( beta.*cosh(beta*L)+(alpha-p11).*sinh(beta*L) );

%G11h =   ( exp(alpha*l).* ( beta.*cosh(beta*(l-L)) + (alpha-p22).*sinh(beta*(l-L))) ) ./ ...
%         ( beta .* cosh(beta*L)-(alpha-p22).*sinh(beta*L) ); dzia�a tylko dla skalarnego l
G11h = repmat(C11.',1,length(l)) .* (exp(alpha.'*(l)).*cosh(beta.'*(l-L))) + ...
       repmat(D11.',1,length(l)) .* (exp(alpha.'*(l)).*sinh(beta.'*(l-L))); 

%G12h =   ( p12*exp(alpha*(l-L)).*sinh(beta*l) ) ./ ... 
%         ( beta.*cosh(beta*L) + (alpha-p11).*sinh(beta*L) ); dzia�a tylko dla skalarnego l
G12h = repmat(C12.',1,length(l)) .* (exp(alpha.'*(l-L)).*sinh(beta.'*l));
          
%G21h =   ( p21*exp(alpha.*l).*sinh(beta.*(l-L)) ) ./ ... 
%         ( beta.*cosh(beta*L) - (alpha-p22).*sinh(beta*L) ); dzia�a tylko dla skalarnego l
G21h = repmat(C21.',1,length(l)) .* (exp(alpha.'*(l)).*sinh(beta.'*(l-L)));

%G22h =   ( exp(alpha*(l-L)).* ( beta.*cosh(beta*l) + (alpha-p11).*sinh(beta*l)) ) ./ ...
%         ( beta .* cosh(beta*L)+(alpha-p11).*sinh(beta*L) ); dzia�a tylko dla skalarnego l
G22h = repmat(C22.',1,length(l)) .* (exp(alpha.'*(l-L)).*cosh(beta.'*l)) + ...
       repmat(D22.',1,length(l)) .* (exp(alpha.'*(l-L)).*sinh(beta.'*l)); 

