% ===========================================================================
%
% Skrypt rysujπcy przestrzenne charakterystyki czasowe rurociπgu 
% w oparciu o model opisany roz≥oøonymi transmitancjami operatorowymi
%
% ===========================================================================

close all
clear

Lp=1000;
omega = logspace(-4,2,100000);
t = [0:0.1:60];


% rozk≥ad ciúnieÒ w stanie ustalonym bez wycieku
l1 = 750;
P = 500000;
[dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek]=oblicz_ustalony(P,Lp,l1,0)


G=trans_widm(omega,Lp,q1);


% odpowiedü: skokowa
hG = zeros(1,length(t));

% obliczanie wartoúci wyraøenia podca≥kowego dla kolejnych chwil czasowych
for i = 1:length(t)
   i
   wart_hG = P_sin_omega(G,omega,t(i));
   
   % ca≥kowanie po omega  
   hG(i) = trapz(omega,wart_hG);
end

hG = 2/pi*hG;

% wykres
figure
plot(t(1:301),hG(1:301,101),'k')
xlabel('t(s)')
ylabel('p(t,Lp)')
%title('\fontname{Arial CE}Odpowiedü skokowa rurociπgu')


% bez tarcia
tt=[0 1.66 1.66 3.33 3.33 5 5 6.66 6.66 8.33 8.33 10 10 11.66 11.66 13.33 13.33 15 15 16.66 16.66 18.33 18.33 20 20 21.66 ...
      21.66 23.33 23.33 25 25 26.66 26.66 28.33 28.33 30];
xx=[-pp -pp pp pp -pp -pp pp pp -pp -pp pp pp -pp -pp pp pp -pp -pp pp pp -pp -pp pp pp -pp -pp pp pp -pp -pp ... 
      pp pp -pp -pp pp pp] ;
   hold on 
   plot(tt,xx,'k--')
   
   
      
