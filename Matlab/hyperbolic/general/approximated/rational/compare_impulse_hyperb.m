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
set(0,'DefaultLineLineWidth',1.5)
set(0,'defaultAxesFontSize',14)

init_hyperb

T = 100;
delta_t=0.01;
t = [0:delta_t:T];  % vector of time instants
l = L;              % spatial position for which the response is evaluated
                              
% calculation of the impulse responses for the distributed parameter model
[g11i,g12i,g21i,g22i] = calc_hyperb_distr_impulse(LAMBDA,K,L,t,l,10);
g11i=real(g11i);
g12i=real(g12i);
g21i=real(g21i);
g22i=real(g22i);

g11i(isnan(g11i))=0;
g12i(isnan(g12i))=0;
g21i(isnan(g21i))=0;
g22i(isnan(g22i))=0;

% calculation of transfer function matrix of the approximation model
% [GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
% calculation of the impulse responses of the approximation model
% [g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);

% frequency responses for the approximated models, N=1, N=10, N=100, N=1000.
for i=1:4
    
 switch i    
    case 1
        load G1_col
        s = 'b:';
    case 2
        load G10_col
        s = 'r-.';
    case 3
        load G100_col
        s = 'g--';
    case 4
        load G1000_col
        s = 'm:';
 end

% impulse responses for approximation models

figure(1)
plot(t,g11Ni,s)
hold on

figure(2)
plot(t,g12Ni,s)
hold on

figure(3)
plot(t,g21Ni,s)
hold on

figure(4)
plot(t,g22Ni,s)
hold on

% step responses for approximation models

h11Ns=delta_t*nancumsum(g11Ni);
h12Ns=delta_t*nancumsum(g12Ni);
h21Ns=delta_t*nancumsum(g21Ni);
h22Ns=delta_t*nancumsum(g22Ni);

figure(5)
plot(t,h11Ns,s)
hold on

figure(6)
plot(t,h12Ns,s)
hold on

figure(7)
plot(t,h21Ns,s)
hold on

figure(8)
plot(t,h22Ns,s)
hold on

end

% impulse responses for distributed parameter models


figure(1)
plot(t,g11i,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',14)
ylabel('g_{11}(t)','FontSize',14)
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')
axis([0 20 -0.05 0.65])
grid on

figure(2)
plot(t,g12i,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',14)
ylabel('g_{12}(t)','FontSize',14)
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')
grid on

figure(3)
plot(t,g21i,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',14)
ylabel('g_{21}(t)','FontSize',14)
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')
grid on

figure(4)
plot(t,g22i,'k-')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',14)
ylabel('g_{22}(t)','FontSize',14)
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')
grid on

% step responses for distributed parameter models

h11s=delta_t*nancumsum(g11i);
h12s=delta_t*nancumsum(g12i);
h21s=delta_t*nancumsum(g21i);
h22s=delta_t*nancumsum(g22i);

figure(5)
plot(t,h11s,'k-')

figure(6)
plot(t,h12s,'k-')

figure(7)
plot(t,h21s,'k-')

figure(8)
plot(t,h22s,'k-')
