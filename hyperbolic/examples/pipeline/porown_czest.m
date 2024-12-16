% ===========================================================================
%
% Skrypt porównuj¹cy charakterystyki czêstotliwoœciowe modu³u ruroci¹gu sprawnego i uszkodzonego
% dla transmitancji Gvl(s,Lp) = P(s,Lp)/Q(s,Lp) 
% 
close all
clear all
clc

Lp=1000;	  % d³ugoœæ ruroci¹gu
P = 5E5;      % ciœnienie na wlocie
l1 = 750;  	
% oblicza spadki ciœnieñ, strumienie oraz rezystancje hydrauliczne w stanie ustalonym
[dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,Lp,l1,0);

omega = 0:0.01:16;
G=trans_widm(omega,Lp,q1);   % transmitacja widmowa ruroci¹gu bez uszkodzeñ


% charakterystyki czestotliwoœciowe modu³u dla 3 ró¿nych œrednic wycieku
figure(1)
plot(omega,log10(abs(G)),'k')	% wykres modu³u transmitacji ruroci¹gu bez uszkodzeñ
%plot(omega,unwrap(angle(G)),'k')	% wykres modu³u transmitacji ruroci¹gu bez uszkodzeñ

l1 = 750;
Dr = [0.01 0.02 0.1];    % œrednice wycieku
col = {'k:' 'k-.' 'k--'};

for i=1:3
   Dl=Dr(i); 
   % oblicza spadki ciœnieñ, strumienie oraz rezystancje hydrauliczne w stanie ustalonym
   [dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,Lp,l1,Dl);
   % transmitacja widmowa ruroci¹gu z wyciekiem
   Gl=trans_widm_uszk(omega,l1,Dl,q1,q2,ql);
   figure(1)
   hold on
   plot(omega,log10(abs(Gl)),col{i}) % wykres modu³u transmitacji ruroci¹gu z wyciekiem
   %plot(omega,unwrap(angle(Gl)),col{i}) % wykres modu³u transmitacji ruroci¹gu z wyciekiem
end
xlabel('\omega [rad/s]')
ylabel('|G_v(j\omega)|')

% charakterystyki czestotliwoœciowe modu³u dla 3 ró¿nych lokalizacji wycieku
figure(2)
plot(omega,log10(abs(G)),'k') % wykres modu³u transmitacji ruroci¹gu bez uszkodzeñ
%plot(omega,unwrap(angle(G)),'k') % wykres modu³u transmitacji ruroci¹gu bez uszkodzeñ

Dl = 0.02;   % Dl=D/10
lr = [250 500 750];

for i=1:3
   l1=lr(i);  
   % oblicza spadki ciœnieñ, strumienie oraz rezystancje hydrauliczne w stanie ustalonym
   [dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,Lp,l1,Dl);
   % transmitacja widmowa ruroci¹gu z wyciekiem
   Gl=trans_widm_uszk(omega,l1,Dl,q1,q2,ql);
   figure(2)
   hold on
   plot(omega,log10(abs(Gl)),col{i}) % wykres modu³u transmitacji ruroci¹gu z wyciekiem
   %plot(omega,unwrap(angle(Gl)),col{i}) % wykres modu³u transmitacji ruroci¹gu z wyciekiem
end
xlabel('\omega [rad/s]')
ylabel('|G_v(j\omega)|')

