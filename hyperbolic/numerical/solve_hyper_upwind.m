function [X1,X2] = solve_hyper_upwind(LAM,K,L,T,XI,XB,M,N)

% Funkcja realizuj¹ca numeryczne rozwi¹zanie metod¹ "upwind" uk³adu dwóch równañ 
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
% zaœ dla ujemnej wartoœci lam - dla punktu l=L+Dl.
%
% Funkcja zwraca macierze X1 oraz X2 o rozmiarach (M+2)*(N+1), zawieraj¹c¹ 
% rozwi¹zania uk³adu w punktach wêz³owych, odpowiednio dla funkcji x1(l,t)
% oraz x2(l,t).

X = zeros(2*(M+2),N+1);   % solution matrix for x1(m,n) and x2(m,n),
                          % for m=0,1,...,M+1 and n=0,1,...,N 
                          % m=0 and m=M+1 represent boundary conditions
                          % n=0 represent initial conditions
                          
%X1 = zeros(M+1,N+1);   % macierz rozwi¹zañ dla zmiennej x1(l,t)
%X2 = zeros(M+1,N+1);   % macierz rozwi¹zañ dla zmiennej x2(l,t)

Dl = L/M;     % wartoœæ kroku dyskretnego w dziedzinie zmiennej przestrzennej
Dt = T/N;     % wartoœæ kroku dyskretnego w dziedzinie czasu

lam1 = LAM(1,1);
lam2 = LAM(2,2);

k11 = K(1,1); k12 = K(1,2);
k21 = K(2,1); k22 = K(2,2);

% Courant numbers
c1 = lam1*Dt/Dl;    % the first equation
c2 = lam2*Dt/Dl;    % the second equation

% checking stability condition
if (abs(c1)>1 | abs(c2)>1)    
    disp('Warning! scheme unstable, Courant numbers:')
    disp('c1='); disp(c1)
    disp('c2='); disp(c2)
    return
end    

% initial conditions
X(1:2:end,1) = XI(:,1);  % x1(l,0)
X(2:2:end,1) = XI(:,2);  % x2(l,0)

% boundary conditions
if (lam1>0)         % x1(0,t)
   X(1,:) = XB(1,:);  
elseif (lam1<0)     % x1(L,t)
   X(end-1,:) = XB(1,:);  
end    
 
if (lam2>0)         
   X(2,:) = XB(2,:);    % x2(0,t) 
elseif (lam2<0)
   X(end,:) = XB(2,:);  % x2(L,t)
end    

% numerical solution is calculated based on the following equation
% X(N+1) = U*X(N)
%
% where Uc is a matrix containing 2*M rows and 2*(M+2) columns

% auxiliary matrix UU for constructing the matrix U
UU = zeros(2*M,5);

UU(1:2:end,3) = -abs(c1)+k11*Dt+1;
UU(1:2:end,4) = k12*Dt;
UU(2:2:end,2) = k21*Dt;
UU(2:2:end,3) = -abs(c2)+k22*Dt+1;

if (lam1>0)
    UU(1:2:end,1) = abs(c1);
elseif (lam1<0)     
    UU(1:2:end,5) = abs(c1);
end    
    
if (lam2>0)
    UU(2:2:end,1) = abs(c2);
elseif (lam2<0)     
    UU(2:2:end,5) = abs(c2);
end    

U = spdiags(UU,0:4,2*M,2*(M+2));

% solution
for n=2:N+1
  X(3:end-2,n) = U*X(:,n-1);
end

% boundary extrapolation
if (lam1>0)         
   X(end-1,1:N+1) = X(end-3,1:N+1);  % x1(M+1,n)=x1(M,n)
elseif (lam1<0)     
   X(1,1:N+1) = X(3,1:N+1);    % x1(0,n)=x1(1,n) 
end    
 
if (lam2>0)         
   X(end,1:N+1) = X(end-2,1:N+1); % x2(M+1,n)=x2(M,n)
elseif (lam2<0)
   X(2,1:N+1) = X(4,1:N+1);  % x2(0,n)=x2(1,n)
end    

X1 = X(1:2:end,:);  % solutions for x1(l,t)
X2 = X(2:2:end,:);  % solutions for x2(l,t)
