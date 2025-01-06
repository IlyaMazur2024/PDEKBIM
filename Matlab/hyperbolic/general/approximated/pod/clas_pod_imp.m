
% Aproksymacja uk�adu o parametrach roz�o�onych typu hiperbolicznego 
% z zastosowaniem analizy g��wnych sk�adowych - odpowied� impulsowa

close all
clear 

K=15;                % liczba wektor�w w�asnych brana pod uwag�
load g_con           % odpowied� impulsowa 

delta_t=0.05;
T = 10;
t = [0:delta_t:T];              % vector of time instants
L = 5;
delta_l=0.1;
l = [0:delta_l:L];              % vector of spatial positions

figure(1)           % wykres przestrzenny odpowiedzi impulsowej
[tt,ll] = meshgrid(l,t);
mesh(tt,ll,real(g21a))
%set(gca,'FontSize',14)
%view(71.5,48)
colormap(gray)
xlabel('l [m]')
ylabel('t [s]')
zlabel('y(l,t)')


Y  = g21a';   % macierz odpowiedzi impulsowej G12, 
                        % warunki brzegowe zgodne
                    
[M,N] = size(Y);    % M - liczba wierszy (pozycji przestrzennych l)  
                    % N - liczba kolumn (pr�bek czasowych t)

yx = mean(Y')';        % �rednia czasowa dla pozycji przestrzennej x 

Yx = repmat(yx,1,N);  % macierz �rednich pozycji przestrzennych

V = Y - Yx;         % macierz odpowiedzi skokowej po usuni�ciu �rednich
                    % dla poszczeg�lnych pozycji przestrzennych
                    % (�rednia dla ka�dego wiersza = 0)    

C = 1/(N-1)*(V*V');  % macierz kowariancji                   
C1 = cov(V);         % macierz kowariancji


[KSI,LAMBDA] = eig(C);  % macierz wektor�w w�asnych i diagonalna macierz warto�ci w�asnych 
                        % macierzy kowariancji

KSI_R = KSI(:,M-K+1:M);  % macierz cech, z�o�ona z n wektor�w w�asnych, 
                           % odpowiadaj�cych najwi�kszym warto�ciom w�asnym    

PSI_R = KSI_R'*V;       % dane, sprowadzone do osi zwi�zanych 
                             % z wybranymi wektorami w�asnymi

figure(2)
%plot(KSI_R','k')
plot(l,KSI_R(:,1),'k-',l,KSI_R(:,2),'k--',l,KSI_R(:,3),'k:' )
legend('\xi_1','\xi_2','\xi_3')
xlabel('l [m]')
ylabel('\xi_{ m}')
set(gca,'FontSize',12)

figure(3)
%plot(PSI_R','k')
plot(t,PSI_R(1,:),'k-',t,PSI_R(2,:),'k--',t,PSI_R(3,:),'k:' )
legend('\psi_1','\psi_2','\psi_3')
xlabel('t [s]')
ylabel('\psi_n')
set(gca,'FontSize',12)

% Kompresja - oryginalne dane w postaci macierzy Y o rozmiarze M x N (51 x 101),
% po usuni�ciu �redniej czasowej zapisane w macierzy V o rozmiarze M x N (51 x 101),
% a nast�pnie skompresowane do macierzy PSI_R o rozmiarze K x N (np. 2 x 101)
% zawieraj�cej dane, sprowadzone do n g��wnych osi zwi�zanych z wektorami 
% w�asnymi macierzy kowariancji, odpowiadaj�cym jej warto�ciom w�asnym o
% najwi�kszej warto�ci.
%
% W celu odtworzenia danych nale�y zapami�ta�: 
% - zredukowan� macierz wektor�w w�asnych KSI_R o rozmiarze M x K, zawieraj�c� wybrane wektory 
%   w�asne macierzy kowariancji (stanowi�ce baz� ortonormaln�)
% - zredukowan� macierz wsp�czynnik�w czasowych PSI_R o rozmiarze K x N
%   (np. 2 x 101)
% - wektor yx o wymiarze 1 x M (np. 51 element�w), zawieraj�cy �rednie
%   czasowe dla poszczeg�lnych pozycji przestrzennych.

% wsp�czynnik "kompresji" danych K:

%               M x N
%  K = -------------------------
%       K x N  +  M x K  +  M

Kc = (M*N)/(K*N+M*K+M)

% procentowa "energia" zwi�zana ze m najwi�kszymi warto�ciami w�asnymi

lam = diag(LAMBDA);
P = sum(lam(M-K+1:M))/sum(lam)

Y_REC = KSI_R*PSI_R + Yx;  % odtworzenie danych
% Y_REC1 = KSI_R*KSI_R'*V + Yx;   % inaczej tak

% ca�kowite odzyskanie oryginalnych danych mo�liwe jedynie dla K=M
% w praktyce - jedynie kilka warto�ci w�asnych r�nych od zera

figure(4)
plot(t,Y(6,:),'k')
hold on
plot(t,Y_REC(6,:),'k--')
plot(t,Y(26,:),'k')
plot(t,Y(46,:),'k')
plot(t,Y_REC(26,:),'k--')
plot(t,Y_REC(46,:),'k--')
xlabel('t [s]')
ylabel('g_{21}(t)')
legend('original response','PCA approximation')
axis([0 10 -0.15 0.8])
text(1,0.7,'l = 0.5')
text(3,0.4,'l = 2.5')
text(6,0.3,'l = 4.5')
set(gca,'FontSize',12)

E = Y - Y_REC;  % b��d rekonstrukcji
figure(5)
mesh(l,t,E')  % wykres b��du rekonstrukcji
set(gca,'FontSize',12)
colormap(gray)
view(68.5,44)
%axis([0 60 0 60 -0.04 0.04])
xlabel('l [m]')
ylabel('t [s]')
zlabel('e(l,t)')

figure(6)
LAM=diag(LAMBDA)
LAM=LAM(end:-1:1);
m=[1:M]
bar(m(1:3),LAM(1:3),'k')
hold on
bar(m(4:end),LAM(4:end),'w')
set(gca,'FontSize',12)
colormap(gray)
%axis([0 60 0 60 -0.04 0.04])
xlabel('m')
ylabel('\sigma_m')
axis([0 51 0 0.22])

% �redniokwadratowy b��d aproksymacji
MSE = 1/(M*N) * sumsqr(E)
% norma Frobeniusa macierzy b��du
no = norm(E,'fro')

E_K = sum(LAM(1:K))/sum(LAM(1:M))