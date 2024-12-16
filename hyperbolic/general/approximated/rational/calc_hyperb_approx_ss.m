
function [SYSN] = calc_hyperb_approx_ss(LAMBDA,K,L,N)

% funkcja oblicza macierze stanu sekcyjnego modelu 
% aproksymacyjnego uk³adu hiperbolicznego 2x2 z³o¿onego z N sekcji

dl = L/N;  % section length 
lambda1 = LAMBDA(1,1);
lambda2 = abs(LAMBDA(2,2));
k11  = K(1,1);
k12  = K(1,2);
k21  = K(2,1);
k22  = K(2,2);

% macierz An
An = [-lambda1/dl+k11 k12; k21 -lambda2/dl+k22];
Bn = [lambda1/dl 0; 0 lambda2/dl];
Cn = [1 0; 0 1];
Dn = [0 0; 0 0];

% macierzowa transmitancja sekcji
SYSn = ss(An,Bn,Cn,Dn);

SYSi = SYSn;
% do³¹czanie kolejnych sekcji modelu

for i=2:N
  SYSi = append(SYSi,SYSn);
end

if lambda1>0 && LAMBDA(2,2)>0 % collocated boundary conditions 

% connection matrix for collocated configuration
connections = [ [3:2+2*(N-1)]' [1:2*(N-1)]' ];
% input vector for collocated configuration
inputs = [1 2];
% output vector for collocated configuration
outputs = [2*N-1 2*N];

% section connection for collocated configuration
SYSN = connect(SYSi,connections,inputs,outputs);

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
SYSN = connect(SYSi,connections,inputs,outputs);

end
