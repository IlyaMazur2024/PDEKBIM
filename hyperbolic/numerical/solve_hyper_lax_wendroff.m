function [X1,X2] = solve_hyper_lax_wendroff(LAM,K,L,T,XI,XB,M,N)

% clear all
% close all
% 
% L=5;    % d³ugoœæ uk³adu
% T=20;   % czas symulacji
% M=50;   % liczba sekcji
% N=200;  % liczba chwil czasowych
% 
% LAM = [1 0; 0 0.5];
% 
% K = [-0.0638  0.0638;
%       0.0359 -0.0359];
% 
% XI = zeros(M+2,2);
% 
% % wymuszenie brzegowe x1(0,t) - funkcja Kroneckera
% XB = [zeros(1,N+1); zeros(1,N+1)];
% XB(1,1) = 1;

% Funkcja realizuj¹ca numeryczne rozwi¹zanie metod¹ Laxa-Wendroffa uk³adu dwóch równañ 
% ró¿niczkowych cz¹stkowych typu hiperbolicznego nastêpuj¹cej postaci:
%
% dx(l,t)/dt + LAM*dx(l,t)/dl = K*x(l,t)
%
% gdzie:   x = [x1(l,t)]
%              [x2(l,t)]
%
% dla:   LAM = [lam1  0   ]   oraz  K = [k11  k12]
%              [ 0   lam2 ]             [k21  k22]
%
% czyli:
% dx1(l,t)/dt + lam1*dx1(l,t)/dl = k11*x1(l,t)+k12*x2(l,t)
% dx2(l,t)/dt + lam2*dx2(l,t)/dl = k21*x1(l,t)+k22*x2(l,t)
%
% dla l z przedzia³u [0,L] oraz t z przedzia³u [0,T]
% odpowiednio z krokiem Dl=L/M oraz Dt=T/N
% 
% dla warunków pocz¹tkowych x1(l,0) i x2(l,0) zawartych 
% w macierzy XI o rozmiarach (M+2)*2 
% oraz warunków brzegowych zawartych w macierzy XB 
% o rozmiarach 2*(N+1), przy czym dla dodatniej wartoœci lam 
% warunki brzegowe powinny byæ okreœlone dla punktu l=0, 
% zaœ dla ujemnej wartoœci lam - dla punktu l=L.
%
% Funkcja zwraca macierze X1 oraz X2 o rozmiarach (M+1)*(N+1), zawieraj¹c¹ 
% rozwi¹zania uk³adu w punktach wêz³owych, odpowiednio dla funkcji x1(l,t)
% oraz x2(l,t).

X = zeros(2*(M+2),N+1);   % solution matrix for x1(m,n) and x2(m,n),
                          % for m=0,1,...,M+1 and n=0,1,...,N 
                          % m=0 and m=M+1 represent boundary conditions
                          % n=0 represent initial conditions
                          
Dl = L/M;     % wartoœæ kroku dyskretnego w dziedzinie zmiennej przestrzennej l
Dt = T/N;     % wartoœæ kroku dyskretnego w dziedzinie czasu t

lam1 = LAM(1,1);
lam2 = LAM(2,2);

k11 = K(1,1); k12 = K(1,2);
k21 = K(2,1); k22 = K(2,2);

c1 = lam1*Dt/Dl;    % Courant number for the first equation
c2 = lam2*Dt/Dl;    % Courant number for the second equation


if (abs(c1)>1 | abs(c2)>1)    % checking stability condition 
    disp('Warning! scheme unstable, Courant numbers:')
    disp('c1='); disp(c1)
    disp('c2='); disp(c2)
    return
end    

X(1:2:end,1) = XI(:,1);  % initial conditions x1(l,0)
X(2:2:end,1) = XI(:,2);  % initial conditions x2(l,0)

if (lam1>0)         % boundary condition x1(0,t)
   X(1,:) = XB(1,:);  
elseif (lam1<0)     % boundary condition x1(L,t)
   X(end-1,:) = XB(1,:);  
end    
 
if (lam2>0)         
   X(2,:) = XB(2,:); % boundary condition x2(0,t) 
elseif (lam2<0)
   X(end,:) = XB(2,:);  % boundary condition x2(L,t)
end    

% numerical solution is calculated based on the following equation
% X(N+1) = P*X(N) + DX(N)*Q*X(N)

% creating band matrix P 
%
% auxiliary matrix LL of diagonal elements
LL = zeros(2*M,6);
LL(1:2:(2*M-1),1) = 1/2*c1*(c1+1);
LL(2:2:2*M,1) = 1/2*c2*(c2+1);
LL(1:2:(2*M-1),2) = 0;
LL(2:2:2*M,2) = k21*Dt;
LL(1:2:(2*M-1),3) = 1+k11*Dt-c1*c1;
LL(2:2:2*M,3) = 1+k22*Dt-c2*c2;
LL(1:2:(2*M-1),4) = k12*Dt;
LL(2:2:2*M,4) = 0;
LL(1:2:(2*M-1),5) = 1/2*c1*(c1-1);
LL(2:2:2*M,5) = 1/2*c2*(c2-1);
P = spdiags(LL,0:5,2*M,2*(M+2));

% band matrix Q of the quadratic part
QQ = [ 1/2*k11*c1*Dt  0  1/2*k11*k11*Dt*Dt  k11*k12*Dt*Dt  -1/2*k11*c1*Dt 0;
       1/2*k12*c1*Dt  0  0  1/2*k12*k12*Dt*Dt  -1/2*k12*c1*Dt  0;
       0  1/2*k21*c2*Dt  1/2*k21*k21*Dt*Dt  k21*k22*Dt*Dt  0  -1/2*k21*c2*Dt;
       0  1/2*k22*c2*Dt  0  1/2*k22*k22*Dt*Dt  0  -1/2*k22*c2*Dt];
Q = zeros(4*M,2*(M+2));
for i=1:M
 Q((4*(i-1)+1):(4*(i-1)+4),(2*(i-1)+1):(2*(i-1)+6)) = QQ;
end 
Q = sparse(Q);

for n=2:N+1     % discrete time instances, n=1...N

  DX = zeros(2*M,4*M);

  for i=1:M
    DX(2*i-1,4*i-3:4*i-2) = X(2*i+1:2*i+2,n-1)';  
    DX(2*i,4*i-1:4*i) = X(2*i+1:2*i+2,n-1)';  
  end
 
  % numerical solution for the sample n+1
  X(3:end-2,n) = P*X(:,n-1) + DX*Q*X(:,n-1);
 
  % artificial boundaries (constant extrapolation)
 
  if (lam1>0)         % artificial boundary x1(M+1,n)=x1(M,n)
    X(end-1,n) = X(end-3,n);  
  elseif (lam1<0)     % artificial boundary x1(0,n)=x1(1,n)
    X(1,n) = X(3,n);  
  end    
 
  if (lam2>0)         % artificial boundary x2(M+1,n)=x2(M,n)
    X(end,n) = X(end-2,n); 
  elseif (lam2<0)
    X(2,n) = X(4,n);  % artificial boundary x2(0,n)=x2(1,n)
  end    
 
end

X1 = X(1:2:end,:);  % solutions for x1(l,t)
X2 = X(2:2:end,:);  % solutions for x2(l,t)
