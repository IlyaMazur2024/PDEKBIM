
function [GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N)

% funkcja oblicza macierz transmitancji operatorowych sekcyjnego modelu 
% aproksymacyjnego uk³adu hiperbolicznego 2x2 z³o¿onego z N sekcji

dl = L/N;  % section length 
lambda1 = LAMBDA(1,1);
lambda2 = abs(LAMBDA(2,2));
k11  = K(1,1);
k12  = K(1,2);
k21  = K(2,1);
k22  = K(2,2);

% wspó³czynniki wielomianu w mianownikach transmitancji sekcji 
a1 = lambda1/dl + lambda2/dl - k11 - k22;
a0 = (k11 - lambda1/dl)*(k22 - lambda2/dl)-k21*k12;

% wspó³czynniki wielomianów w licznikach transmitancji sekcji 
b11_1 = lambda1/dl;      b11_0 = lambda1/dl*(lambda2/dl-k22);
b12_0 = k12*lambda2/dl;  b21_0 = k21*lambda1/dl; 
b22_1 = lambda2/dl;      b22_0 = lambda2/dl*(lambda1/dl-k11);

% wektor wspó³czynników w mianowniku transmitancji
m   = [1 a1 a0];
% wektory wspó³czynników w liczniku transmitancji
l11 = [b11_1 b11_0];
l12 = [b12_0];
l21 = [b21_0];
l22 = [b22_1 b22_0];

% transmitancje operatorowe poszczególnych torów sekcji
G11n = tf(l11,m);
G12n = tf(l12,m);
G21n = tf(l21,m);
G22n = tf(l22,m);

% macierzowa transmitancja sekcji
Gn = [G11n G12n; G21n G22n];

if N==1
    GN=Gn;
    return;
end    

Gi = Gn;
% do³¹czanie kolejnych sekcji modelu

for i=2:N
  Gi = append(Gi,Gn);
end

if lambda1>0 && LAMBDA(2,2)>0 % collocated boundary conditions 

% connection matrix for collocated configuration
connections = [ [3:2+2*(N-1)]' [1:2*(N-1)]' ];
% input vector for collocated configuration
inputs = [1 2];
% output vector for collocated configuration
outputs = [2*N-1 2*N];

% section connection for collocated configuration
GN = connect(Gi,connections,inputs,outputs);

elseif lambda1>0 && LAMBDA(2,2)<0 % anti-collocated boundary conditions

% connection matrix for anti-collocated configuration
connections = zeros(2*(N-1),2);
connections(1:2:2*N-3,:) = [ [3:2:2*N-1]' [1:2:2*N-3]' ];
connections(2:2:2*N-2,:) = [ [2:2:2*N-2]' [4:2:2*N]' ];

% input vector for anti-collocated configuration
inputs = [1 2*N];

% output vector for anti-collocated configuration
outputs = [2*N-1 2];
%outputs = [1:2*N]; 

% section connection for anti-collocated configuration
GN = connect(Gi,connections,inputs,outputs);

end





