
% Aproksymacja uk³adu o parametrach roz³o¿onych typu hiperbolicznego 
% z zastosowaniem analizy g³ównych sk³adowych - odpowiedŸ impulsowa

close all
clear 

K=15;                % liczba wektorów w³asnych brana pod uwagê
load g_con           % odpowiedŸ impulsowa 

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
                    % N - liczba kolumn (próbek czasowych t)

yx = mean(Y')';        % œrednia czasowa dla pozycji przestrzennej x 

Yx = repmat(yx,1,N);  % macierz œrednich pozycji przestrzennych

V = Y - Yx;         % macierz odpowiedzi skokowej po usuniêciu œrednich
                    % dla poszczególnych pozycji przestrzennych
                    % (œrednia dla ka¿dego wiersza = 0)    

C = 1/(N-1)*(V*V');  % macierz kowariancji                   
C1 = cov(V);         % macierz kowariancji


[KSI,LAMBDA] = eig(C);  % macierz wektorów w³asnych i diagonalna macierz wartoœci w³asnych 
                        % macierzy kowariancji

KSI_R = KSI(:,M-K+1:M);  % macierz cech, z³o¿ona z n wektorów w³asnych, 
                           % odpowiadaj¹cych najwiêkszym wartoœciom w³asnym    

PSI_R = KSI_R'*V;       % dane, sprowadzone do osi zwi¹zanych 
                             % z wybranymi wektorami w³asnymi

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
% po usuniêciu œredniej czasowej zapisane w macierzy V o rozmiarze M x N (51 x 101),
% a nastêpnie skompresowane do macierzy PSI_R o rozmiarze K x N (np. 2 x 101)
% zawieraj¹cej dane, sprowadzone do n g³ównych osi zwi¹zanych z wektorami 
% w³asnymi macierzy kowariancji, odpowiadaj¹cym jej wartoœciom w³asnym o
% najwiêkszej wartoœci.
%
% W celu odtworzenia danych nale¿y zapamiêtaæ: 
% - zredukowan¹ macierz wektorów w³asnych KSI_R o rozmiarze M x K, zawieraj¹c¹ wybrane wektory 
%   w³asne macierzy kowariancji (stanowi¹ce bazê ortonormaln¹)
% - zredukowan¹ macierz wspó³czynników czasowych PSI_R o rozmiarze K x N
%   (np. 2 x 101)
% - wektor yx o wymiarze 1 x M (np. 51 elementów), zawieraj¹cy œrednie
%   czasowe dla poszczególnych pozycji przestrzennych.

% wspó³czynnik "kompresji" danych K:

%               M x N
%  K = -------------------------
%       K x N  +  M x K  +  M

Kc = (M*N)/(K*N+M*K+M)

% procentowa "energia" zwi¹zana ze m najwiêkszymi wartoœciami w³asnymi

lam = diag(LAMBDA);
P = sum(lam(M-K+1:M))/sum(lam)

Y_REC = KSI_R*PSI_R + Yx;  % odtworzenie danych
% Y_REC1 = KSI_R*KSI_R'*V + Yx;   % inaczej tak

% ca³kowite odzyskanie oryginalnych danych mo¿liwe jedynie dla K=M
% w praktyce - jedynie kilka wartoœci w³asnych ró¿nych od zera

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

E = Y - Y_REC;  % b³¹d rekonstrukcji
figure(5)
mesh(l,t,E')  % wykres b³êdu rekonstrukcji
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

% œredniokwadratowy b³¹d aproksymacji
MSE = 1/(M*N) * sumsqr(E)
% norma Frobeniusa macierzy b³êdu
no = norm(E,'fro')

E_K = sum(LAM(1:K))/sum(LAM(1:M))