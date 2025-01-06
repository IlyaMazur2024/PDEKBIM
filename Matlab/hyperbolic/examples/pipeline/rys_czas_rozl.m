% ===========================================================================
%
% Skrypt rysuj�cy przestrzenne charakterystyki czasowe ruroci�gu 
% w oparciu o model opisany roz�o�onymi transmitancjami operatorowymi
%
% ===========================================================================

close all
clear

Lp=1000;
omega = logspace(-4,2,100000);
t = [0:0.1:60];
x = [0:10:1000];

% rozk�ad ci�nie� w stanie ustalonym bez wycieku
l1 = 750;
P = 500000;
[dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek]=oblicz_ustalony(P,Lp,l1,0)

G=trans_widm(omega,x,q1);

% odpowied�: skokowa
hG = zeros(length(t),length(x));


% obliczanie warto�ci wyra�enia podca�kowego dla kolejnych chwil czasowych i kolejnych punkt�w
  for i = 1:length(t)
    i
    wart_hG = P_sin_omega(G,omega,t(i));
    % ca�kowanie po omega  
    hG(i,:) = trapz(omega,wart_hG',1);
  end

hG = 2/pi*hG;

% wykres
figure
t = [0:0.1:30];
[xx,yy] = meshgrid(x,t);
zz = hG(1:301,:);
surf(xx,yy,zz)
view(82.5,54)
colormap(gray)
shading interp
view(96,44)
xlabel('l [m]')
ylabel('t [s]')
zlabel('p(t,l)')
