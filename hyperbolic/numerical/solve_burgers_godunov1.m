function [X1,X2] = solve_burgers_godunov1(L,T,XI,M,N)

% Funkcja realizuj¹ca numeryczne rozwi¹zanie metod¹ Godunowa uk³adu 
% równañ ró¿niczkowych cz¹stkowych Burgersa:
%
% dx1(l,t)/dt + x1(l,t)*dx1(l,t)/dl = 0
% dx2(l,t)/dt + x2(l,t)*dx2(l,t)/dl = 0
%
% dla l z przedzia³u [0,L] oraz t z przedzia³u [0,T]
% odpowiednio z krokiem Dl=L/M oraz Dt=T/N
% 
% dla warunku pocz¹tkowego x1(l,0) i x2(l,0) zawartych 
% w macierzy XI o rozmiarach (M+2)*2 
%
% Funkcja zwraca macierze X1 oraz X2 o rozmiarach (M+2)*(N+1), zawieraj¹c¹ 
% rozwi¹zania uk³adu w punktach wêz³owych, odpowiednio dla funkcji x1(l,t)
% oraz x2(l,t).

X = zeros(2*(M+2),N+1);   % solution matrix for x1(m,n) and x2(m,n),
                          % for m=0,1,...,M+1 and n=0,1,...,N 
                          % m=0 and m=M+1 represent boundary conditions
                          % n=0 represent initial conditions
                          
Dl = L/M;     % wartoœæ kroku dyskretnego w dziedzinie zmiennej przestrzennej
Dt = T/N;     % wartoœæ kroku dyskretnego w dziedzinie czasu

% initial conditions
X(1:2:end,1) = XI(:,1);  % x1(l,0)
X(2:2:end,1) = XI(:,2);  % x2(l,0)


% numerical solution is calculated based on the following equation
% X(n+1) = U*X^2(n) + X(n)
%
% where U is a matrix containing 2*M rows and 2*(M+2) columns

% auxiliary matrix UU for constructing the matrix U
UU = zeros(2*M,3);

UU(1:end,1) = Dt/(2*Dl);
UU(1:end,3) = -Dt/(2*Dl);

U = spdiags(UU,0:2,2*M,2*(M+2));

% solution
for n=2:N+1
  X(3:end-2,n) = U*X(:,n-1).^2+X(3:end-2,n-1);
end

% % boundary extrapolation
% if (lam1>0)         
%    X(end-1,1:N+1) = X(end-3,1:N+1);  % x1(M+1,n)=x1(M,n)
% elseif (lam1<0)     
%    X(1,1:N+1) = X(3,1:N+1);    % x1(0,n)=x1(1,n) 
% end    
%  
% if (lam2>0)         
%    X(end,1:N+1) = X(end-2,1:N+1); % x2(M+1,n)=x2(M,n)
% elseif (lam2<0)
%    X(2,1:N+1) = X(4,1:N+1);  % x2(0,n)=x2(1,n)
% end    

X1 = X(1:2:end,:);  % solutions for x1(l,t)
X2 = X(2:2:end,:);  % solutions for x2(l,t)
