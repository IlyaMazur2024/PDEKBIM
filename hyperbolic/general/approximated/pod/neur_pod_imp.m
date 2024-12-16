
% Aproksymacja uk³adu o parametrach roz³o¿onych typu hiperbolicznego 
% z zastosowaniem liniowej sieci neuronowej - odpowiedŸ impulsowa

close all
clear 

N=15;                % liczba wektorów w³asnych brana pod uwagê

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
set(gca,'FontSize',12)
%view(71.5,48)
colormap(gray)
xlabel('l [m]')
ylabel('t [s]')
zlabel('y(l,t)')

P = g21a';
T = g21a';

[c,r] = size(P);   % c - liczba pozycji przestrzennych


net = newff(minmax(P),[N c],{'purelin','purelin'})
net.trainFcn='trainlm';
net.trainParam.epochs = 40;
net.trainParam.goal = 1E-6;
net = train(net,P,T)

% pierwsza warstwa sieci (realizuj¹ca kompresjê)
net1=newlin(minmax(P),N,0,0.001);
net1.IW{1}=net.IW{1};
net1.b{1}=net.b{1};

v = sim(net1,P);          % dane czasowe sprowadzone do N wymiarów

W = net.IW{1}              % wagi pierwszej warstwy
figure(2)
%plot(l,W(1,:),'k-',l,W(2,:),'k--',l,W(3,:),'k:' )
xlabel('l [m]')
ylabel('w_i')
legend('w_1','w_2','w_3')
set(gca,'FontSize',12)

figure(3)
%plot(t,v(1,:),'k-',t,v(2,:),'k--',t,v(3,:),'k:' )
xlabel('t [s]')
ylabel('v_i')
legend('v_1','v_2','v_3')
set(gca,'FontSize',12)

% druga warstwa sieci (realizuj¹ca dekompresjê)
net2=newlin(minmax(v),51,0,0.001);     
net2.IW{1}=net.LW{2};
net2.b{1}=net.b{2};

y = sim(net2,v)         % odtworzenie (dekompresja) danych

% Kompresja - oryginalne dane przestrzenno-czasowe w postaci macierzy P 
% o rozmiarze c x r (51 x 101) 
% zosta³y skompresowane do macierzy odpowiedzi pierwszej warstwy y1 
% o wymiarach N x r (5 x 101) zawieraj¹cej dane, sprowadzone do N g³ównych osi.

% Dodatkowo w celu umo¿liwienia rekonstrukcji danych nale¿y zapamiêtaæ: 
% macierz wspó³czynników wagowych drugiej warstwy sieci o rozmiarach c x N (np. 51 x 5)
% oraz wektor wspó³czynników progowych tej warstwy o wymiarze c x 1 (51 elementów).

% wspó³czynnik kompresji K:

%               r x c
%  K = -------------------------
%       N x r  +  c x N  +  c

K = (r*c)/(N*r+c*N+c)

% ca³kowite odzyskanie oryginalnych danych mo¿liwe jedynie dla N=c

figure(4)               % comparative plot for l=0.5, l=2.5, l=4.5
plot(t,P(6,:),'k')
hold on
plot(t,y(6,:),'k--')
plot(t,P(26,:),'k')
plot(t,P(46,:),'k')
plot(t,y(26,:),'k--')
plot(t,y(46,:),'k--')
xlabel('t [s]')
ylabel('g_{21}(t)')
legend('original response','FF-PCA approximation')
axis([0 10 -0.15 0.8])
text(1,0.7,'l = 0.5')
text(3,0.4,'l = 2.5')
text(6,0.3,'l = 4.5')
set(gca,'FontSize',12)

E = P - y;  % b³¹d rekonstrukcji
figure(5)
mesh(l,t,E')  % wykres b³êdu rekonstrukcji
set(gca,'FontSize',12)
colormap(gray)
view(68.5,44)
%axis([0 60 0 60 -0.04 0.04])
xlabel('l [m]')
ylabel('t [s]')
zlabel('e(l,t)')

[M,N] = size(P);
% œredniokwadratowy b³¹d aproksymacji
SSE = 1/(M*N) * sumsqr(E)
% norma Frobeniusa macierzy b³êdu
no = norm(E,'fro')

% xlabel('iterations')
% ylabel('MSAE')

