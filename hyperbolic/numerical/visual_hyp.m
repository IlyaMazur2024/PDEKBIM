
clear all
close all

L=5;    % d³ugoœæ uk³adu
T=20;   % czas symulacji
M=50;   % liczba sekcji
N=400;  % liczba chwil czasowych
lam1 = 1;
lam2 = 0.2;
k11 =  0.05; k12 = 0.05;
k21 =  0.05; k22 = 0.05;

LAM = [lam1 0; 0 lam2];
%K = [-0.0638  0.0638;
%      0.0359 -0.0359];
K = [k11 k12;
     k21 k22];

s1=(k22*lam1-k11*lam2+2*sqrt(-k12*k21*lam1*lam2))/(lam1-lam2)
s2=(k22*lam1-k11*lam2-2*sqrt(-k12*k21*lam1*lam2))/(lam1-lam2)
     
% L=1;    % d³ugoœæ uk³adu
% T=75;   % czas symulacji
% M=100;   % liczba sekcji
% N=9375;  % liczba chwil czasowych
% LAM = [1 0; 0 0.5];
% K = [0 0;
%      0 0];

% warunki pocz¹tkowe - zerowe 
XI = zeros(M+2,2);
%XI(:,1)=sin(0.1*(1:M+2));
%XI(5:10,1)=1;

% wymuszenie brzegowe x1(0,t) - funkcja Kroneckera
XB = [zeros(1,N+1); zeros(1,N+1)];
XB(1,1) = 1;

% wymuszenie brzegowe x1(0,t) - funkcja Heaviside'a
%XB = [ones(1,N+1); zeros(1,N+1)];

% wymuszenie brzegowe x1(0,t) - impuls prostok¹tny
%XB = [zeros(1,N+1); zeros(1,N+1)];
%XB(1,1:N/4) = 1;

% wymuszenie brzegowe x1(0,t) - fala sinusoidalna
%XB = [sin(0.05*(1:N+1)); zeros(1,N+1)];

[X1,X2] = solve_hyper_upwind(LAM,K,L,T,XI,XB,M,N);
%[X1,X2] = solve_hyper_lax_wendroff(LAM,K,L,T,XI,XB,M,N);
%[X1,X2,c1,c2] = solve_hyper_leapfrog(LAM,K,L,T,XI,XB,M,N);
%[X1,X2] = solve_burgers_godunov1(L,T,XI,M,N);

Dl = L/M;     
Dt = T/N; 

% wizualizacja rozwi¹zania
[TT,LL]=meshgrid([0:Dt:T],[0:Dl:L+Dl]);
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

