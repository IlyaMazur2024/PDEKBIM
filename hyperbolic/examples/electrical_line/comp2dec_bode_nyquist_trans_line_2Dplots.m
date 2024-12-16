% ===========================================================================
%
% Skrypt porównuj¹cy charakterystyki czêstotliwoœciowe 2D Bodego i Nyquista 
% uk³adu hiperbolicznego z jedn¹ zmienn¹ przestrzenn¹ w oparciu o model 
% opisany transmitancjami operatorowymi dla uk³adu o parametrach
% roz³o¿onych, warunki brzegowe zgodne i przeciwstawne
%

close all
clear all

trans_line_init                % inicjalizacja wartoœci modelu wymiennika
omega = logspace(-1.5,1.5,10000); % wektor pulsacji dla wymiennika


% wyznaczenie odpowiedzi czêstotliwoœciowych uk³adu oryginalnego (sprzê¿onego)

l = Lp;                        % charakterystyki dla punktu koñcowego
[G11wc,G12wc,G21wc,G22wc,G11wi,G12wi,dum,dum] = freq_resp_trans_line(omega,l,R,L,G,C);

l = 0;                         % charakterystyki dla punktu pocz¹tkowego
[dum,dum,dum,dum,dum,dum,G21wi,G22wi] = freq_resp_trans_line(omega,l,R,L,G,C);

LAMBDA = diag([-1/sqrt(L*C), 1/sqrt(L*C)]);
K = [ -R/(2*L)-G/(2*C),   G/(2*C)-R/(2*L);
       G/(2*C)-R/(2*L),  -R/(2*L)-G/(2*C)];


% wyznaczenie odpowiedzi czêstotliwoœciowych uk³adu po rozsprzê¿eniu (s³abo sprzê¿onego)

l = Lp;                        % charakterystyki dla punktu koñcowego
[G11c,G12c,G21c,G22c,G11i,G12i,dum,dum] = freq_resp_decoupled(omega,l,LAMBDA,K);

l = 0;                         % charakterystyki dla punktu pocz¹tkowego
[dum,dum,dum,dum,dum,dum,G21i,G22i] = freq_resp_decoupled(omega,l,LAMBDA,K);  
   

% ===========================================================================
% Warunki brzegowe zgodne

% Charakterystyki Bodego

figure(1)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11wc)),'k')
hold on
plot(log10(omega),20*log10(abs(G11c)),'k:')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_1_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G11wc)),'k')
hold on
plot(log10(omega),unwrap(angle(G11c)),'k:')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_1(j\omega)')
%axis([-4 1 -60 10])


figure(2)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12wc)),'k')
hold on
plot(log10(omega),20*log10(abs(G12c)),'k:')
xlabel('lg(\omega)')
ylabel('M(L,\omega) [dB]')
%title('Amplitude plot for G_1_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12wc)),'k')
hold on
plot(log10(omega),unwrap(angle(G12c)),'k:')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_2(j\omega)')
%axis([-4 1 -60 10])
legend('G_{w12}(L,\omega)','G_{12}(L,\omega)')

figure(3)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21wc)),'k')
hold on
plot(log10(omega),20*log10(abs(G21c)),'k:')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_2_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21wc)),'k')
hold on
plot(log10(omega),unwrap(angle(G21c)),'k:')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_1(j\omega)')
%axis([-4 1 -60 10])

figure(4)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22wc)),'k')
hold on
plot(log10(omega),20*log10(abs(G22c)),'k:')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_2_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22wc)),'k')
hold on
plot(log10(omega),unwrap(angle(G22c)),'k:')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_2(j\omega)')
%axis([-4 1 -60 10])



% Charakterystyki amplitudowo-fazowe Nyquista
%
figure(5)
plot(real(G11wc),imag(G11wc),'k');
hold on
plot(real(G11c),imag(G11c),'k:');
xlabel('Re\langleG_1_1(L,j\omega)\rangle')
ylabel('Im\langleG_1_1(L,j\omega)\rangle')
%title('Nyquist plot for G_1_1(j\omega)')

figure(6)
plot(real(G12wc),imag(G12wc),'k');
hold on
plot(real(G12c),imag(G12c),'k:');
xlabel('Re\langleG_1_2(L,j\omega)\rangle')
ylabel('Im\langleG_1_2(L,j\omega)\rangle')
%title('Nyquist plot for G_1_2(j\omega)')

figure(7)
plot(real(G21wc),imag(G21wc),'k');
hold on
plot(real(G21c),imag(G21c),'k:');
xlabel('Re\langleG_2_1(L,j\omega)\rangle')
ylabel('Im\langleG_2_1(L,j\omega)\rangle')
%title('Nyquist plot for G_2_1(j\omega)')

figure(8)
plot(real(G22wc),imag(G22wc),'k');
hold on
plot(real(G22c),imag(G22c),'k:');
xlabel('Re\langleG_2_2(L,j\omega)\rangle')
ylabel('Im\langleG_2_2(L,j\omega)\rangle')
%title('Nyquist plot for G_2_2(j\omega)')

pause
close all

% ===========================================================================
% Warunki brzegowe przeciwstawne

% Charakterystyki Bodego

figure(1)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11wi)),'k')
hold on
plot(log10(omega),20*log10(abs(G11i)),'k:')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
title('Amplitude plot for G_1_1(j\omega)')
subplot(2,1,2)
%plot(log10(omega),unwrap(angle(G11wi)),'k')
hold on
plot(log10(omega),unwrap(angle(G11i)),'k:')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_1(j\omega)')
%axis([-4 legend('G_{w12}(L_p,\omega)','G_{12}(L_p,\omega)')

figure(2)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12wi)),'k')
hold on
plot(log10(omega),20*log10(abs(G12i)),'k:')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_1_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12wi)),'k')
hold on
plot(log10(omega),unwrap(angle(G12i)),'k:')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_2(j\omega)')
%axis([-4 1 -60 10])
legend('G_{w12}(L,\omega)','G_{12}(L,\omega)')

figure(3)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21wi)),'k')
hold on
plot(log10(omega),20*log10(abs(G21i)),'k:')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_2_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21wi)),'k')
hold on
plot(log10(omega),unwrap(angle(G21i)),'k:')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_1(j\omega)')
%axis([-4 1 -60 10])

figure(4)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22wi)),'k')
hold on
plot(log10(omega),20*log10(abs(G22i)),'k:')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_2_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22wi)),'k')
hold on
plot(log10(omega),unwrap(angle(G22i)),'k:')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_2(j\omega)')
%axis([-4 1 -60 10])

% Charakterystyki amplitudowo-fazowe Nyquista
%
figure(5)
plot(real(G11wi),imag(G11wi),'k');
hold on
plot(real(G11i),imag(G11i),'k:');
xlabel('Re\langleG_1_1(L,j\omega)\rangle')
ylabel('Im\langleG_1_1(L,j\omega)\rangle')
%title('Nyquist plot for G_1_1(j\omega)')

figure(6)
%plot(real(G12wi(1:6500)),imag(G12wi(1:6500)),'k');
hold on
plot(real(G12i),imag(G12i),'k:');
xlabel('Re\langleG_1_2(L,j\omega)\rangle')
ylabel('Im\langleG_1_2(L,j\omega)\rangle')
%title('Nyquist plot for G_1_2(j\omega)')

figure(7)
plot(real(G21wi),imag(G21wi),'k');
hold on
plot(real(G21i),imag(G21i),'k:');
xlabel('Re\langleG_2_1(L,j\omega)\rangle')
ylabel('Im\langleG_2_1(L,j\omega)\rangle')
%title('Nyquist plot for G_2_1(j\omega)')

figure(8)
plot(real(G22wi),imag(G22wi),'k');
hold on
plot(real(G22i),imag(G22i),'k:');
xlabel('Re\langleG_2_2(L,j\omega)\rangle')
ylabel('Im\langleG_2_2(L,j\omega)\rangle')
%title('Nyquist plot for G_2_2(j\omega)')

