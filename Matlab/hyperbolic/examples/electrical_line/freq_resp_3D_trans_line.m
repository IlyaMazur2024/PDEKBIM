% ===========================================================================
%
% Skrypt rysuj¹cy charakterystyki czêstotliwoœciowe 3D 
% uk³adu hiperbolicznego z jedn¹ zmienn¹ przestrzenn¹ w oparciu o model 
% opisany transmitancjami operatorowymi dla uk³adu o parametrach
% roz³o¿onych, warunki pocz¹tkowe przeciwne
clc
close all
clear all

trans_line_init                   % inicjalizacja wartoœci modelu 
omega = logspace(-1.5,1.5,100); % wektor pulsacji 

l = [0:10:Lp];                  % wektor zmiennych przestrzennych

[G11wc,G12wc,G21wc,G22wc,G11wi,G12wi,G21wi,G22wi] = freq_resp_trans_line(omega,l,R,L,G,C);

% Charakterystyki przestrzenne (3D) czêœci rzeczywistych i urojonych
% transmitancji widmowych - warunki brzegow przeciwstawne
%
set(0,'DefaultAxesFontSize', 12)

figure(1)
[xx,yy] = meshgrid(l,log10(omega));
zz = real(G11wi);
mesh(xx,yy,zz)
view(78,62)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('Re\langleG_w_1_1(l,j\omega)\rangle')
%title('\fontname{Arial CE}Czêœæ rzeczywista transmitancji G_1_1(l,j\omega) dla przeciwpr¹du')

figure(2)
[xx,yy] = meshgrid(l,log10(omega));
zz = imag(G11wi);
mesh(xx,yy,zz)
view(78,62)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('Im\langleG_w_1_1(l,j\omega)')
%title('\fontname{Arial CE}Czêœæ urojona transmitancji G_1_1(l,j\omega) dla przeciwpr¹du')

figure(3)
[xx,yy] = meshgrid(l,log10(omega));
zz = real(G12wi);
mesh(xx,yy,zz)
view(72,54)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('Re\langleG_w_1_2(l,j\omega)\rangle')
%title('\fontname{Arial CE}Czêœæ rzeczywista transmitancji G_1_2(l,j\omega) dla przeciwpr¹du')

figure(4)
[xx,yy] = meshgrid(l,log10(omega));
zz = imag(G12wi);
mesh(xx,yy,zz)
view(75,40)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('Im\langleG_w_1_2(l,j\omega)\rangle')
%title('\fontname{Arial CE}Czêœæ urojona transmitancji G_1_2(l,j\omega) dla przeciwpr¹du')

figure(5)
[xx,yy] = meshgrid(l,log10(omega));
zz = real(G21wi);
mesh(xx,yy,zz)
view(72,54)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('Re\langleG_w_2_1(l,j\omega)\rangle')
%title('\fontname{Arial CE}Czêœæ rzeczywista transmitancji G_2_1(l,j\omega) dla przeciwpr¹du')

figure(6)
[xx,yy] = meshgrid(l,log10(omega));
zz = imag(G21wi);
mesh(xx,yy,zz)
view(75,40)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('Im\langleG_w_2_1(l,j\omega)\rangle')
%title('\fontname{Arial CE}Czêœæ urojona transmitancji G_2_1(l,j\omega) dla przeciwpr¹du')

figure(7)
[xx,yy] = meshgrid(l,log10(omega));
zz = real(G22wi);
mesh(xx,yy,zz)
view(78,62)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('Re\langleG_w_2_2(l,j\omega)\rangle')
%title('\fontname{Arial CE}Czêœæ rzeczywista transmitancji G_2_2(l,j\omega) dla przeciwpr¹du')

figure(8)
[xx,yy] = meshgrid(l,log10(omega));
zz = imag(G22wi);
mesh(xx,yy,zz)
view(78,62)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('Im\langleG_w_2_2(l,j\omega)\rangle')
%title('\fontname{Arial CE}Czêœæ urojona transmitancji G_2_2(l,j\omega) dla przeciwpr¹du')


% =========================================================================

% Przestrzenne charakterystyki czêstotliwoœciowe Bodego 
%

figure(9)
[xx,yy] = meshgrid(l,log10(omega));
zz = 20*log10(abs(G11wi));
mesh(xx,yy,zz)
view(58,22)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('A_v(l,\omega) [dB]')
%title('\fontname{Arial CE}Charakterystyka amplitudowa G_1_1(l,j\omega) dla przeciwpr¹du')

figure(10)
[xx,yy] = meshgrid(l,log10(omega));
zz = unwrap(angle(G11wi));
mesh(xx,yy,zz)
view(60,52)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('\phi(l,\omega [rad]')
%title('\fontname{Arial CE}Charakterystyka fazowa G_1_1(l,j\omega) dla przeciwpr¹du')


figure(11)
[xx,yy] = meshgrid(l,log10(omega));
zz = 20*log10(abs(G12wi));
mesh(xx,yy,zz)
view(77,40)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('A_v(l,\omega) [dB]')
%title('\fontname{Arial CE}Charakterystyka amplitudowa G_1_2(l,j\omega) dla przeciwpr¹du')

figure(12)
[xx,yy] = meshgrid(l,log10(omega));
zz = unwrap(angle(G12wi));
mesh(xx,yy,zz)
view(76,56)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('\phi(l,\omega) [rad]')
%title('\fontname{Arial CE}Charakterystyka fazowa G_1_2(l,j\omega) dla przeciwpr¹du')

figure(13)
[xx,yy] = meshgrid(l,log10(omega));
zz = 20*log10(abs(G21wi));
mesh(xx,yy,zz)
view(77,40)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('A_v(l,\omega) [dB]')
%title('\fontname{Arial CE}Charakterystyka amplitudowa G_2_1(l,j\omega) dla przeciwpr¹du')

figure(14)
[xx,yy] = meshgrid(l,log10(omega));
zz = unwrap(angle(G21wi));
mesh(xx,yy,zz)
view(76,56)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('\phi(l,\omega) [rad]')
%title('\fontname{Arial CE}Charakterystyka fazowa G_2_1(l,j\omega) dla przeciwpr¹du')

figure(15)
[xx,yy] = meshgrid(l,log10(omega));
zz = 20*log10(abs(G22wi));
mesh(xx,yy,zz)
view(63,60)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('A_v(l,\omega) [dB]')
%title('\fontname{Arial CE}Charakterystyka amplitudowa G_2_2(l,j\omega) dla przeciwpr¹du')

figure(16)
[xx,yy] = meshgrid(l,log10(omega));
zz = unwrap(angle(G22wi));
mesh(xx,yy,zz)
view(77,46)
colormap(cool)
xlabel('l')
ylabel('lg(\omega)')
zlabel('\phi(l,\omega) [rad]')
%title('\fontname{Arial CE}Charakterystyka fazowa G_2_2(l,j\omega) dla przeciwpr¹du')


