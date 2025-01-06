% =========================================================================
%
% m-file name: compare_impulse_2D_responses.m
%
% Purpose: comparison of the impulse responses of the 2 x 2 hyperblic 
% system calculated based on the analytical expressions and 
%
% =========================================================================

clc
clear all 
close all

trans_line_init  

delta_t=0.01;
T = 30;
t = [0:delta_t:T];              % vector of time instants
l = 1000;                        % current spatial position

[g11a,g12a,g21a,g22a]=numerical_impulse_responses_line(t,l,R,L,G,C);

% while(any(isnan(g11a)))
%     g11a(isnan(g11a)) = g11a(find(isnan(g11a))-1);
% end
% while(any(isnan(g12a)))
%     g12a(isnan(g12a)) = g12a(find(isnan(g12a))-1);
% end
% while(any(isnan(g21a)))
%     g21a(isnan(g21a)) = g21a(find(isnan(g21a))-1);
% end
% while(any(isnan(g22a)))
%     g22a(isnan(g22a)) = g22a(find(isnan(g22a))-1);
% end

% g11a(isnan(g11a))=0;
% g12a(isnan(g12a))=0;
% g21a(isnan(g21a))=0;
% g22a(isnan(g22a))=0;

% step responses

figure(1)
%plot(t,g11n,'b')
hold on
plot(t,g11a,'k:')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('g_{11}(l,t)','FontSize',16)
%legend('numerical impulse response','analytical impulse response')
grid on

figure(2)
%plot(t,g12n,'b')
hold on
plot(t,g12a,'k:')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('g_{12}(l,t)','FontSize',16)
%legend('numerical impulse response','analytical impulse response')
grid on

figure(3)
%plot(t,g21n,'b')
hold on
plot(t,g21a,'k:')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('g_{21}(l,t)','FontSize',16)
%legend('numerical impulse response','analytical impulse response')
grid on

figure(4)
%plot(t,g22n,'b')
hold on
plot(t,g22a,'k:')
set(gca,'FontSize',14)
xlabel('t [s]','FontSize',16)
ylabel('g_{22}(l,t)','FontSize',16)
%legend('numerical impulse response','analytical impulse response')
grid on

% step responses

% h11n=delta_t*nancumsum(g11n);
 h11a=delta_t*nancumsum(g11a);
% h12n=delta_t*nancumsum(g12n);
 h12a=delta_t*nancumsum(g12a);
% h21n=delta_t*nancumsum(g21n);
 h21a=delta_t*nancumsum(g21a);
% h22n=delta_t*nancumsum(g22n);
 h22a=delta_t*nancumsum(g22a);

 figure(5)
% plot(t,h11n,'b')
% hold on
 plot(t,h11a,'r')
% legend('numerical step response','analytical step response')
% 
 figure(6)
% plot(t,h12n,'b')
% hold on
 plot(t,h12a,'r')
% legend('numerical step response','analytical step response')
% 
 figure(7)
% plot(t,h21n,'b')
% hold on
 plot(t,h21a,'r')
% legend('numerical step response','analytical step response')
% 
 figure(8)
% plot(t,h22n,'b')
% hold on
 plot(t,h22a,'r')
% legend('numerical step response','analytical step response')
