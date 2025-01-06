% ===========================================================================
%
% Skrypt rysuj¹cy charakterystyki czêstotliwoœciowe ruroci¹gu dla l=Lr
% dla transmitancji G(s,L) = P(s,L)/Q(s,L) = - Z*tanh(gamma*l)
% 

close all
clear

Lp=1000;
%omega = 0:0.01:200;
omega = 0:0.01:20;

% rozk³ad ciœnieñ w stanie ustalonym bez wycieku
l1 = 750;
P = 500000;
[dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek]=oblicz_ustalony(P,Lp,l1,0)


G=trans_widm(omega,Lp,q1);


% Charakterystyki czêstotliwoœciowe Bodego 
%
figure
%plot(log10(omega),20*log10(abs(G)),'k')
plot(omega,20*log10(abs(G)),'k')
xlabel('\omega [rad/s]')
ylabel('log|(G_v(j\omega))|')
%title('\fontname{Arial CE}Charakterystyka amplitudowa G_v(j\omega) ruroci¹gu')
%axis([0 20 4 8])

figure
%plot(log10(omega),unwrap(angle(G)),'k')
plot(omega,unwrap(angle(G)),'k')
xlabel('\omega [rad/s]')
ylabel('\phi(\omega) [rad]')
%title('\fontname{Arial CE}Charakterystyka fazowa G_v(j\omega) ruroci¹gu')
%axis([0 20 -5 -1])


