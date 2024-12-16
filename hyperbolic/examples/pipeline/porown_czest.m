% ===========================================================================
%
% Skrypt por�wnuj�cy charakterystyki cz�stotliwo�ciowe modu�u ruroci�gu sprawnego i uszkodzonego
% dla transmitancji Gvl(s,Lp) = P(s,Lp)/Q(s,Lp) 
% 
close all
clear all
clc

Lp=1000;	  % d�ugo�� ruroci�gu
P = 5E5;      % ci�nienie na wlocie
l1 = 750;  	
% oblicza spadki ci�nie�, strumienie oraz rezystancje hydrauliczne w stanie ustalonym
[dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,Lp,l1,0);

omega = 0:0.01:16;
G=trans_widm(omega,Lp,q1);   % transmitacja widmowa ruroci�gu bez uszkodze�


% charakterystyki czestotliwo�ciowe modu�u dla 3 r�nych �rednic wycieku
figure(1)
plot(omega,log10(abs(G)),'k')	% wykres modu�u transmitacji ruroci�gu bez uszkodze�
%plot(omega,unwrap(angle(G)),'k')	% wykres modu�u transmitacji ruroci�gu bez uszkodze�

l1 = 750;
Dr = [0.01 0.02 0.1];    % �rednice wycieku
col = {'k:' 'k-.' 'k--'};

for i=1:3
   Dl=Dr(i); 
   % oblicza spadki ci�nie�, strumienie oraz rezystancje hydrauliczne w stanie ustalonym
   [dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,Lp,l1,Dl);
   % transmitacja widmowa ruroci�gu z wyciekiem
   Gl=trans_widm_uszk(omega,l1,Dl,q1,q2,ql);
   figure(1)
   hold on
   plot(omega,log10(abs(Gl)),col{i}) % wykres modu�u transmitacji ruroci�gu z wyciekiem
   %plot(omega,unwrap(angle(Gl)),col{i}) % wykres modu�u transmitacji ruroci�gu z wyciekiem
end
xlabel('\omega [rad/s]')
ylabel('|G_v(j\omega)|')

% charakterystyki czestotliwo�ciowe modu�u dla 3 r�nych lokalizacji wycieku
figure(2)
plot(omega,log10(abs(G)),'k') % wykres modu�u transmitacji ruroci�gu bez uszkodze�
%plot(omega,unwrap(angle(G)),'k') % wykres modu�u transmitacji ruroci�gu bez uszkodze�

Dl = 0.02;   % Dl=D/10
lr = [250 500 750];

for i=1:3
   l1=lr(i);  
   % oblicza spadki ci�nie�, strumienie oraz rezystancje hydrauliczne w stanie ustalonym
   [dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,Lp,l1,Dl);
   % transmitacja widmowa ruroci�gu z wyciekiem
   Gl=trans_widm_uszk(omega,l1,Dl,q1,q2,ql);
   figure(2)
   hold on
   plot(omega,log10(abs(Gl)),col{i}) % wykres modu�u transmitacji ruroci�gu z wyciekiem
   %plot(omega,unwrap(angle(Gl)),col{i}) % wykres modu�u transmitacji ruroci�gu z wyciekiem
end
xlabel('\omega [rad/s]')
ylabel('|G_v(j\omega)|')

