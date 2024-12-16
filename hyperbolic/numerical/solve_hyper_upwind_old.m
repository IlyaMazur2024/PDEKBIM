% function [X1,X2] = solve_hyper_upwind_old(LAM,K,L,T,XI,XB,M,N)

L=5;    % d³ugoœæ uk³adu
T=20;   % czas symulacji
M=50;   % liczba sekcji
N=200;  % liczba chwil czasowych

LAM = [0.5 0; 0 0.5];

%K = [-0.0638  0.0638;
%      0.0359 -0.0359];

K = [0 0;
     0 0];

% warunki pocz¹tkowe - zerowe 
XI = zeros(M+1,2);

% wymuszenie brzegowe x1(0,t) - funkcja Kroneckera
XB = [zeros(1,N+1); zeros(1,N+1)];
XB(1,1) = 1;

% wymuszenie brzegowe x1(0,t) - funkcja Heaviside'a
%XB = [ones(1,N+1); zeros(1,N+1)];

% wymuszenie brzegowe x1(0,t) - impuls prostok¹tny
XB = [zeros(1,N+1); zeros(1,N+1)];
XB(1,1:N/4) = 1;

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
% w macierzy XI o rozmiarach (M+1)*2 
% oraz warunków brzegowych zawartych w macierzy XB 
% o rozmiarach 2*(N+1), przy czym dla dodatniej wartoœci lam 
% warunki brzegowe powinny byæ okreœlone dla punktu l=0, 
% zaœ dla ujemnej wartoœci lam - dla punktu l=L.
%
% Funkcja zwraca macierze X1 oraz X2 o rozmiarach (M+1)*(N+1), zawieraj¹c¹ 
% rozwi¹zania uk³adu w punktach wêz³owych, odpowiednio dla funkcji x1(l,t)
% oraz x2(l,t).

X1 = zeros(M+1,N+1);   % macierz rozwi¹zañ dla zmiennej x1(l,t)
X2 = zeros(M+1,N+1);   % macierz rozwi¹zañ dla zmiennej x2(l,t)

Dl = L/M;     % wartoœæ kroku dyskretnego w dziedzinie zmiennej przestrzennej
Dt = T/N;     % wartoœæ kroku dyskretnego w dziedzinie czasu

lam1 = LAM(1,1);
lam2 = LAM(2,2);

k11 = K(1,1); k12 = K(1,2);
k21 = K(2,1); k22 = K(2,2);

c1 = lam1*Dt/Dl;    % Courant number for the first equation
c2 = lam2*Dt/Dl;    % Courant number for the second equation

if (abs(c1)>1 | abs(c2)>1)    % stability condition
    disp('Warning! scheme unstable, Courant numbers:')
    disp('c1='); disp(c1)
    disp('c2='); disp(c2)
end    

X1(:,1) = XI(:,1);  % initial conditions x1(l,0)
X2(:,1) = XI(:,2);  % initial conditions x2(l,0)

if (lam1>0)         % boundary conditions x1(0,t)
    X1(1,:) = XB(1,:);  
elseif (lam1<0)
    X1(M+1,:) = XB(1,:);  
end    

if (lam2>0)         % boundary conditions x2(0,t)
    X2(1,:) = XB(2,:);  
elseif (lam2<0)
    X2(M+1,:) = XB(2,:);  
end    

A1 = abs(c1)*diag(ones(M+1,1)) + (-abs(c1)+k11*Dt+1)*diag(ones(M,1),1) ;
A1 = A1(1:M,:);
B1 = k12*Dt*diag(ones(M,1));

A2 = abs(c2)*diag(ones(M+1,1)) + (-abs(c2)+k22*Dt+1)*diag(ones(M,1),1);
A2 = A2(1:M,:);
B2 = k21*Dt*diag(ones(M,1));

if (lam1>0)
  for n=2:N+1           % discrete time instances
     X1(2:M+1,n) = A1*X1(1:M+1,n-1) + B1*X2(2:M+1,n-1);
  end
elseif (lam1<0)  
  for n=2:N+1           % discrete time instances
     X1(M:-1:1,n) = A1*X1(M+1:-1:1,n-1) + B1*X2(M:-1:1,n-1);
  end  
end

if (lam2>0)
  for n=2:N+1           % discrete time instances
     X2(2:M+1,n) = A2*X2(1:M+1,n-1) + B2*X1(2:M+1,n-1);
  end
elseif (lam2<0)  
  for n=2:N+1           % discrete time instances
     X2(M:-1:1,n) = A2*X2(M+1:-1:1,n-1) + B2*X1(M:-1:1,n-1);
  end  
end

[TT,LL]=meshgrid([0:Dt:T],[0:Dl:L]);
figure(1)
mesh(TT,LL,X1)
xlabel('t')
ylabel('l')
title('x_1(l,t)')

figure(2)
mesh(TT,LL,X2)
xlabel('t')
ylabel('l')
title('x_2(l,t)')

