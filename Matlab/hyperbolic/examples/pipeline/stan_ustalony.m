
% Skrypt wyznaczaj¹cy rozk³ad ciœnienia oraz przep³ywu wzd³u¿ osi ruroci¹gu
% 
%
clear all
close all
clc

L = 1000;
l1 = 750;
ll = [0 L];
ll1 = [0 l1];
ll2 = [l1 L];
P = 500000;

% rozk³ad ciœnieñ w stanie ustalonym bez wycieku

[dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek]=oblicz_ustalony(P,L,l1,0)

p1 = [P P-dp1];
p2 = [P-dp1 P-dp1-dp2];
qq1 = [q1 q1];
qq2 = [q2 q2];

figure(1)
plot(ll1,p1,'k--')
hold on
plot(ll2,p2,'k--')
xlabel('l [m]')
ylabel('p [Pa]')
axis([0 1000 460000 520000])
plot([750 750],[460000 520000],'b:')

figure(2)
plot(ll1,qq1,'k--')
hold on
plot(ll2,qq2,'k--')
xlabel('l [m]')
ylabel('q [kg/s]')
plot([750 750],[5 35],'b:')



% rozk³ad ciœnieñ w stanie ustalonym z wyciekiem
Dl=0.02;		% œrednica wycieku
[dp1u,dp2u,dplu,dpku,q1u,q2u,qlu,R1u,R2u,Rlu,Rku,Re1u,Re2u,Relu,Reku]=oblicz_ustalony(P,L,l1,Dl)

p1 = [P P-dp1u];
p2 = [P-dp1u P-dp1u-dp2u];
qq1 = [q1u q1u];
qq2 = [q2u q2u];
qql = [qlu qlu];


figure(1)
plot(ll1,p1,'r')
hold on
plot(ll2,p2,'r')

figure(2)
plot(ll1,qq1,'r')
hold on
plot(ll2,qq2,'r')
plot(ll,qql,'r')


