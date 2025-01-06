function [X1,X2,c1,c2] = solve_hyper_leapfrog(LAM,K,L,T,XI,XB,M,N)

% Funkcja realizuj¹ca numeryczne rozwi¹zanie metod¹ leap-frog uk³adu dwóch równañ 
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

% two-step scheme - requires a one-level scheme to get started
[X1,X2] = solve_hyper_upwind(LAM,K,L,Dt,XI,XB(:,1:2),M,1);
X(1:2:end,1:2) = X1;  
X(2:2:end,1:2) = X2;  

% numerical solution is calculated based on the following 
% matrix difference equation:
% X(n+1) = P*X(n) + X(n-1)

% creating band matrix P 
%
% auxiliary matrix PP of diagonal elements

PP = zeros(2*M,5);

PP(1:2:end,1) = c1;
PP(2:2:end,1) = c2;
PP(2:2:end,2) = 2*k21*Dt;
PP(1:2:end,3) = 2*k11*Dt;
PP(2:2:end,3) = 2*k22*Dt;
PP(1:2:end,4) = 2*k12*Dt;
PP(1:2:end,5) = -c1;
PP(2:2:end,5) = -c2;

P = spdiags(PP,0:4,2*M,2*(M+2));

for n=3:N+1     % discrete time instances, n=1...N

 % numerical solution for the sample n+1
  X(3:end-2,n) = P*X(:,n-1) + X(3:end-2,n-2);
 
  % artificial boundaries (constant extrapolation)
 
   if (lam1>0)         % artificial boundary x1(M+1,n)=x1(M,n-1)
     X(end-1,n) = X(end-3,n-1);  
   elseif (lam1<0)     % artificial boundary x1(0,n)=x1(1,n-1)
     X(1,n) = X(3,n-1);  
   end    
  
   if (lam2>0)         % artificial boundary x2(M+1,n)=x2(M,n-1)
     X(end,n) = X(end-2,n-1); 
   elseif (lam2<0)
     X(2,n) = X(4,n-1);  % artificial boundary x2(0,n)=x2(1,n-1)
   end    
 
end

X1 = X(1:2:end,:);  % solutions for x1(l,t)
X2 = X(2:2:end,:);  % solutions for x2(l,t)
