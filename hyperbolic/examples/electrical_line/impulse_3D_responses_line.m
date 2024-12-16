% =========================================================================
%
% m-file name: compare_impulse-3D_responses.m
%
% Purpose: comparison of the spatiotemporal impulse responses of the 
% 2 x 2 hyperblic system calculated based on the analytical expressions 
% and numerical approach
% 
% =========================================================================

clear all 
close all

trans_line_init  

delta_t=0.1;
T = 30;
t = [0:delta_t:T];              % vector of time instants

L = 1000;
delta_l=10;
l = [0:delta_l:L];                  % vector of spatial positions

[g11a,g12a,g21a,g22a]= numerical_impulse_responses_line(t,l,R,L,G,C);

g11a=real(g11a);
g12a=real(g12a);
g21a=real(g21a);
g22a=real(g22a);

g11a(isnan(g11a))=0;
g12a(isnan(g12a))=0;
g21a(isnan(g21a))=0;
g22a(isnan(g22a))=0;

figure(1)
[ll,tt] = meshgrid(l,t);
mesh(ll,tt,real(g11a))
colormap(cool)
set(gca,'FontSize',16)
xlabel('l [m]','FontSize',16)
ylabel('t [s]','FontSize',16)
zlabel('g_{11}(l,t)','FontSize',16)

figure(2)
[ll,tt] = meshgrid(l,t);
mesh(ll,tt,real(g12a))
colormap(cool)
set(gca,'FontSize',14)
xlabel('l [m]','FontSize',16)
ylabel('t [s]','FontSize',16)
zlabel('g_{12}(l,t)','FontSize',16)


figure(3)
[ll,tt] = meshgrid(l,t);
mesh(ll,tt,real(g21a))
colormap(cool)
set(gca,'FontSize',14)
xlabel('l [m]','FontSize',16)
ylabel('t [s]','FontSize',16)
zlabel('g_{21}(l,t)','FontSize',16)


figure(4)
[ll,tt] = meshgrid(l,t);
mesh(ll,tt,real(g22a))
colormap(cool)
set(gca,'FontSize',14)
xlabel('l [m]','FontSize',16)
ylabel('t [s]','FontSize',16)
zlabel('g_{22}(l,t)','FontSize',16)

