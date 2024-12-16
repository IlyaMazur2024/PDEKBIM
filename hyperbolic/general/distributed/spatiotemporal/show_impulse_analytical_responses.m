% =========================================================================
%
% m-file name: 
%
% Purpose: comparison of the impulse responses of the 2 x 2 hyperblic 
% system calculated based on the analytical expressions and   
%
% =========================================================================

clc
clear all 
close all

init_hyperb

T = 500;
delta_t=0.01;
t = [0:delta_t:T];  % vector of time instants
l = L;              % spatial position for which the response is evaluated
                              
% calculation of the impulse responses for the distributed parameter model
[g11i,g12i,g21i,g22i] = analytical_impulse_responses(LAMBDA,K,L,t,l,10);
g11i=real(g11i);
g12i=real(g12i);
g21i=real(g21i);
g22i=real(g22i);

g11i(isnan(g11i))=0;
g12i(isnan(g12i))=0;
g21i(isnan(g21i))=0;
g22i(isnan(g22i))=0;

% impulse responses for distributed parameter models


figure(1)
plot(t,g11i,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('g_{11}(t)','FontSize',16)
grid on

figure(2)
plot(t,g12i,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('g_{12}(t)','FontSize',16)
grid on

figure(3)
plot(t,g21i,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('g_{21}(t)','FontSize',16)
grid on

figure(4)
plot(t,g22i,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('g_{22}(t)','FontSize',16)
grid on

% step responses for distributed parameter models

h11s=delta_t*nancumsum(g11i);
h12s=delta_t*nancumsum(g12i);
h21s=delta_t*nancumsum(g21i);
h22s=delta_t*nancumsum(g22i);

figure(5)
plot(t,h11s,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('h_{11}(t)','FontSize',16)
grid on

figure(6)
plot(t,h12s,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('h_{12}(t)','FontSize',16)
grid on

figure(7)
plot(t,h21s,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('h_{21}(t)','FontSize',16)
grid on


figure(8)
plot(t,h22s,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('h_{22}(t)','FontSize',16)
grid on
