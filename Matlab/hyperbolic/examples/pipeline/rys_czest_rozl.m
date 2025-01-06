% ===========================================================================
%
% Skrypt rysuj¹cy przestrzenne charakterystyki czêstotliwoœciowe ruroci¹gu
% dla transmitancji G(s,L) = P(s,L)/Q(s,L) = - Z*tanh(gamma*l)
% 

close all
clear

% rozk³ad ciœnieñ w stanie ustalonym bez wycieku
L = 1000;
l1 = 750;
P = 500000;
[dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek]=oblicz_ustalony(P,L,l1,0)


omega = 0:0.01:20;
x = [0:1:1000];


G=trans_widm(omega,x,q1);


% Przestrzenne charakterystyki czêstotliwoœciowe Bodego 
%
figure(1)
%[xx,yy] = meshgrid(x,log10(omega));
[xx,yy] = meshgrid(x,omega);
zz = 20*log10(abs(G));
mesh(xx,yy,zz)
view(64,22)
colormap(hsv)
shading interp;   
xlabel('l [m]')
ylabel('log(\omega)')
ylabel('\omega [rad/s]')
zlabel('L(\omega,l) [dB]')
%title('\fontname{Arial CE}Charakterystyka amplitudowa G(j\omega,l) ruroci¹gu')

figure(2)
x = [5:1:1000];
G=trans_widm(omega,x,q1);
%[xx,yy] = meshgrid(x,log10(omega));
[xx,yy] = meshgrid(x,omega);
zz = unwrap(angle(G));
mesh(xx,yy,zz)
view(60,52)
colormap(hsv)
shading interp;   
xlabel('l [m]')
%ylabel('log(\omega)')
ylabel('\omega [rad/s]')
zlabel('\phi(\omega,l) [rad]')
%title('\fontname{Arial CE}Charakterystyka fazowa G(j\omega,x) ruroci¹gu')

