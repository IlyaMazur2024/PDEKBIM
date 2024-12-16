
% Aproksymacja uk�adu o parametrach roz�o�onych typu hiperbolicznego 
% z zastosowaniem sieci neuronowej typu PCA uczonej uog�lnion� metod� Hebba
% (odpowied� impulsowa)

close all
clear 
clc

K=15;                % liczba wektor�w w�asnych brana pod uwag�

load g_con           % odpowied� impulsowa 

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
                    % N - liczba kolumn (pr�bek czasowych t)

yx = mean(Y')';        % �rednia czasowa dla pozycji przestrzennej x 

Yx = repmat(yx,1,N);  % macierz �rednich pozycji przestrzennych

Ym = Y - Yx;         % macierz odpowiedzi skokowej po usuni�ciu �rednich
                    % dla poszczeg�lnych pozycji przestrzennych
                    % (�rednia dla ka�dego wiersza = 0)    

%[W, errvals]=gha(Ym,K,'rate',1,'niter',101);   % macierz wsp�czynnik�w wagowych sieci PCA
                        % czyli (kolumnowych) wektor�w w�asnych 
                       
                        
% =========================================================================uczenie sieci PCA

W = 0.1*randn(M,K);  % losowa macierz pocz�tkowych wsp�czynnik�w wagowych
IT = 5000;           % liczba iteracji
eta = 0.01;             % wsp�czynnik pr�dko�ci uczenia
SE = zeros(1,IT+1);

tic
for i = 1:IT
    
    % wyznaczanie b��du 
    
    Y_R = W*W'*Ym;    % 
    E = Ym - Y_R;     % 
    SE(i) = 1/(M*N) * sumsqr(E);   % b��d przed i-t� iteracj�  
    
    for j = 1:N
      x = Ym(:,j);    % kolejny wzorzec ucz�cy - warto�ci zmiennej procesowej
                     % dla wszystkich pozycji przestrzennych w kolejnej 
                     % chwili czasowej
      y = W'*x;      % podanie wzorca ucz�cego na M wej�� sieci
                     % y - wektor K sygna��w wyj�ciowych   
      
      W = eta*x*y' + W - eta*W*tril(y*y');  % modyfikacja wag - algorytm GHA
            
    end
end
toc

Y_R = W*W'*Ym;    % 
E = Ym - Y_R;     % 
SE(i+1) = 1/(M*N) * sumsqr(E);   % b��d po zako�czeniu ostatniej iteracji  


% =========================================================================      

                        
V = W'*Ym;               % macierz zawieraj�ca odpowiedzi sieci na 
                        % kolejne wektory wej�ciowe 

PHI_R = W;              % macierz wag sieci odpowiada zredukowanej macierzy PHI
                        % (ortogonalnych wektor�w zmiennej przestrzennej)
                        
PSI_R = V;              % macierz odpowiedzi sieci odpowiedzi sieci odpowiada 
                        % zredukowanej macierzy PSI (ortogonalnych wektor�w
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
% po usuni�ciu �redniej czasowej zapisane w macierzy V o rozmiarze M x N (51 x 101),
% a nast�pnie skompresowane do macierzy PSI_R o rozmiarze K x N (np. 2 x 101)
% zawieraj�cej dane, sprowadzone do n g��wnych osi zwi�zanych z wektorami 
% w�asnymi macierzy kowariancji, odpowiadaj�cym jej warto�ciom w�asnym o
% najwi�kszej warto�ci.
%
% W celu odtworzenia danych nale�y zapami�ta�: 
% - zredukowan� macierz wektor�w w�asnych PHI_R o rozmiarze M x K, zawieraj�c� wybrane wektory 
%   w�asne macierzy kowariancji (stanowi�ce baz� ortonormaln�)
% - zredukowan� macierz wsp�czynnik�w czasowych PSI_R o rozmiarze K x N
%   (np. 2 x 101)
% - wektor yx o wymiarze 1 x M (np. 51 element�w), zawieraj�cy �rednie
%   czasowe dla poszczeg�lnych pozycji przestrzennych.

% wsp�czynnik "kompresji" danych K:

%               M x N
%  K = -------------------------
%       K x N  +  M x K  +  M

c = (M*N)/(K*N+M*K+M)

% procentowa "energia" zwi�zana ze m najwi�kszymi warto�ciami w�asnymi


Y_REC = PHI_R*PSI_R + Yx;  % odtworzenie danych
% Y_REC1 = PHI_R*PHI_R'*V + Yx;   % inaczej tak

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
legend('original response','GHA-PCA approximation')
axis([0 10 -0.15 0.8])
text(1,0.7,'l = 0.5')
text(3,0.4,'l = 2.5')
text(6,0.3,'l = 4.5')
set(gca,'FontSize',12)

E = Y - Y_REC;  % b��d rekonstrukcji
figure(5)
mesh(tt(1:50,:),ll(1:50,:),E(:,1:50)')  % wykres b��du rekonstrukcji
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

% �redniokwadratowy b��d aproksymacji
SSE = 1/(M*N) * sumsqr(E)
% norma Frobeniusa macierzy b��du
no = norm(E,'fro')
