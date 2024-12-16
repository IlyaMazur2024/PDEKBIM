
% Aproksymacja uk³adu o parametrach roz³o¿onych typu hiperbolicznego 
% z zastosowaniem sieci neuronowej typu PCA uczonej uogólnion¹ metod¹ Hebba
% (odpowiedŸ impulsowa)

close all
clear 
clc

K=15;                % liczba wektorów w³asnych brana pod uwagê

load g_con           % odpowiedŸ impulsowa 

delta_t=0.05;
T = 10;
t = [0:delta_t:T];              % vector of time instants
L = 5;
delta_l=0.1;
l = [0:delta_l:L];              % vector of spatial positions

figure(1)            % wykres przestrzenny odpowiedzi impulsowej
[tt,ll] = meshgrid(l,t);
mesh(tt,ll,real(g21a))
set(gca,'FontSize',12)
%view(71.5,48)
colormap(gray)
xlabel('l')
ylabel('t')
zlabel('y(l,t)')

Y  = g21a';           % macierz odpowiedzi impulsowej G12, 
                      % warunki brzegowe zgodne
                   
[M,N] = size(Y);    % M - liczba wierszy (pozycji przestrzennych l)  
                    % N - liczba kolumn (próbek czasowych t)

yx = mean(Y')';        % œrednia czasowa dla pozycji przestrzennej x 

Yx = repmat(yx,1,N);  % macierz œrednich pozycji przestrzennych

Ym = Y - Yx;         % macierz odpowiedzi skokowej po usuniêciu œrednich
                    % dla poszczególnych pozycji przestrzennych
                    % (œrednia dla ka¿dego wiersza = 0)    

%[W, errvals]=gha(Ym,K,'rate',1,'niter',101);   % macierz wspó³czynników wagowych sieci PCA
                        % czyli (kolumnowych) wektorów w³asnych 
                       
                        
% =========================================================================uczenie sieci PCA

W = 0.1*randn(M,K);  % losowa macierz pocz¹tkowych wspó³czynników wagowych
IT = 5000;           % liczba iteracji
eta = 0.01;             % wspó³czynnik prêdkoœci uczenia
SE = zeros(1,IT+1);

tic
for i = 1:IT
    
    % wyznaczanie b³êdu 
    
    Y_R = W*W'*Ym;    % 
    E = Ym - Y_R;     % 
    SE(i) = 1/(M*N) * sumsqr(E);   % b³¹d przed i-t¹ iteracj¹  
    
    for j = 1:N
      x = Ym(:,j);    % kolejny wzorzec ucz¹cy - wartoœci zmiennej procesowej
                     % dla wszystkich pozycji przestrzennych w kolejnej 
                     % chwili czasowej
      y = W'*x;      % podanie wzorca ucz¹cego na M wejœæ sieci
                     % y - wektor K sygna³ów wyjœciowych   
      
      W = eta*x*y' + W - eta*W*tril(y*y');  % modyfikacja wag - algorytm GHA
            
    end
end
toc

Y_R = W*W'*Ym;    % 
E = Ym - Y_R;     % 
SE(i+1) = 1/(M*N) * sumsqr(E);   % b³¹d po zakoñczeniu ostatniej iteracji  


% =========================================================================      

                        
V = W'*Ym;               % macierz zawieraj¹ca odpowiedzi sieci na 
                        % kolejne wektory wejœciowe 

PHI_R = W;              % macierz wag sieci odpowiada zredukowanej macierzy PHI
                        % (ortogonalnych wektorów zmiennej przestrzennej)
                        
PSI_R = V;              % macierz odpowiedzi sieci odpowiedzi sieci odpowiada 
                        % zredukowanej macierzy PSI (ortogonalnych wektorów
                        % zmiennej czasowej)
                    
                
figure(2)
plot(l,W,'k')
%plot(l,W(:,1),'k-',l,W(:,2),'k--',l,W(:,3),'k:' )
xlabel('l')
ylabel('w_i')
%legend('w_1','w_2','w_3')
set(gca,'FontSize',12)

figure(3)
plot(t,V,'k')
%plot(t,V(1,:),'k-',t,V(2,:),'k--',t,V(3,:),'k:' )
xlabel('t')
ylabel('v_i')
%legend('v_1','v_2','v_3')
set(gca,'FontSize',12)

% Kompresja - oryginalne dane w postaci macierzy Y o rozmiarze M x N (51 x 101),
% po usuniêciu œredniej czasowej zapisane w macierzy V o rozmiarze M x N (51 x 101),
% a nastêpnie skompresowane do macierzy PSI_R o rozmiarze K x N (np. 2 x 101)
% zawieraj¹cej dane, sprowadzone do n g³ównych osi zwi¹zanych z wektorami 
% w³asnymi macierzy kowariancji, odpowiadaj¹cym jej wartoœciom w³asnym o
% najwiêkszej wartoœci.
%
% W celu odtworzenia danych nale¿y zapamiêtaæ: 
% - zredukowan¹ macierz wektorów w³asnych PHI_R o rozmiarze M x K, zawieraj¹c¹ wybrane wektory 
%   w³asne macierzy kowariancji (stanowi¹ce bazê ortonormaln¹)
% - zredukowan¹ macierz wspó³czynników czasowych PSI_R o rozmiarze K x N
%   (np. 2 x 101)
% - wektor yx o wymiarze 1 x M (np. 51 elementów), zawieraj¹cy œrednie
%   czasowe dla poszczególnych pozycji przestrzennych.

% wspó³czynnik "kompresji" danych K:

%               M x N
%  K = -------------------------
%       K x N  +  M x K  +  M

c = (M*N)/(K*N+M*K+M)

% procentowa "energia" zwi¹zana ze m najwiêkszymi wartoœciami w³asnymi


Y_REC = PHI_R*PSI_R + Yx;  % odtworzenie danych
% Y_REC1 = PHI_R*PHI_R'*V + Yx;   % inaczej tak

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
legend('original response','GHA-PCA approximation')
axis([0 10 -0.15 0.8])
text(1,0.7,'l = 0.5')
text(3,0.4,'l = 2.5')
text(6,0.3,'l = 4.5')
set(gca,'FontSize',12)

E = Y - Y_REC;  % b³¹d rekonstrukcji
figure(5)
mesh(tt(1:50,:),ll(1:50,:),E(:,1:50)')  % wykres b³êdu rekonstrukcji
set(gca,'FontSize',12)
colormap(gray)
view(68.5,44)
%axis([0 60 0 60 -0.04 0.04])
xlabel('l')
ylabel('t')
zlabel('e(l,t)')

figure(6)
semilogy(1:IT+1,SE,'k-')
set(gca,'FontSize',12)
xlabel('Iterations')
ylabel('MSAE')
axis([0 5000 0.00005 0.05])

% œredniokwadratowy b³¹d aproksymacji
SSE = 1/(M*N) * sumsqr(E)
% norma Frobeniusa macierzy b³êdu
no = norm(E,'fro')
