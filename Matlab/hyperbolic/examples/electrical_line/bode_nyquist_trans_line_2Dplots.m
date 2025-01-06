% ===========================================================================
%
% Skrypt rysuj¹cy charakterystyki czêstotliwoœciowe 2D Bodego i Nyquista 
% uk³adu hiperbolicznego z jedn¹ zmienn¹ przestrzenn¹ w oparciu o model 
% opisany transmitancjami operatorowymi dla uk³adu o parametrach
% roz³o¿onych, warunki brzegowe zgodne i przeciwstawne

close all
clear all

trans_line_init                % inicjalizacja wartoœci modelu wymiennika
omega = logspace(-1.5,1.5,10000); % wektor pulsacji dla wymiennika

l = Lp;                        % charakterystyki dla punktu koñcowego
[G11c,G12c,G21c,G22c,G11i,G12i,dum,dum] = freq_resp_trans_line(omega,l,R,L,G,C);

l = 0;                        % charakterystyki dla punktu pocz¹tkowego
[dum,dum,dum,dum,dum,dum,G21i,G22i] = freq_resp_trans_line(omega,l,R,L,G,C);

% ===========================================================================
% Warunki brzegowe zgodne

% Charakterystyki Bodego

figure(1)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11c)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_1_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G11c)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_1(j\omega)')
%axis([-4 1 -60 10])

figure(2)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12c)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_1_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12c)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_2(j\omega)')
%axis([-4 1 -60 10])

figure(3)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21c)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_2_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21c)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_1(j\omega)')
%axis([-4 1 -60 10])

figure(4)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22c)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_2_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22c)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_2(j\omega)')
%axis([-4 1 -60 10])



% Charakterystyki amplitudowo-fazowe Nyquista
%
figure(5)
plot(real(G11c),imag(G11c),'r');
xlabel('Re\langleG_1_1(L,j\omega)\rangle')
ylabel('Im\langleG_1_1(L,j\omega)\rangle')
%title('Nyquist plot for G_1_1(j\omega)')

figure(6)
plot(real(G12c),imag(G12c),'r');
xlabel('Re\langleG_1_2(L,j\omega)\rangle')
ylabel('Im\langleG_1_2(L,j\omega)\rangle')
%title('Nyquist plot for G_1_2(j\omega)')

figure(7)
plot(real(G21c),imag(G21c),'r');
xlabel('Re\langleG_2_1(L,j\omega)\rangle')
ylabel('Im\langleG_2_1(L,j\omega)\rangle')
%title('Nyquist plot for G_2_1(j\omega)')

figure(8)
plot(real(G22c),imag(G22c),'r');
xlabel('Re\langleG_2_2(L,j\omega)\rangle')
ylabel('Im\langleG_2_2(L,j\omega)\rangle')
%title('Nyquist plot for G_2_2(j\omega)')

pause

% ===========================================================================
% Warunki brzegowe przeciwstawne

% Charakterystyki Bodego

figure(1)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11i)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_1_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G11i)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_1(j\omega)')
%axis([-4 1 -60 10])

figure(2)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12i)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_1_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12i)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_2(j\omega)')
%axis([-4 1 -60 10])

figure(3)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21i)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_2_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21i)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_1(j\omega)')
%axis([-4 1 -60 10])

figure(4)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22i)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
%title('Amplitude plot for G_2_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22i)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_2(j\omega)')
%axis([-4 1 -60 10])

% Charakterystyki amplitudowo-fazowe Nyquista
%
figure(5)
plot(real(G11i),imag(G11i),'r');
xlabel('Re\langleG_1_1(L,j\omega)\rangle')
ylabel('Im\langleG_1_1(L,j\omega)\rangle')
%title('Nyquist plot for G_1_1(j\omega)')

figure(6)
plot(real(G12i),imag(G12i),'r');
xlabel('Re\langleG_1_2(L,j\omega)\rangle')
ylabel('Im\langleG_1_2(L,j\omega)\rangle')
%title('Nyquist plot for G_1_2(j\omega)')

figure(7)
plot(real(G21i),imag(G21i),'r');
xlabel('Re\langleG_2_1(L,j\omega)\rangle')
ylabel('Im\langleG_2_1(L,j\omega)\rangle')
%title('Nyquist plot for G_2_1(j\omega)')

figure(8)
plot(real(G22i),imag(G22i),'r');
xlabel('Re\langleG_2_2(L,j\omega)\rangle')
ylabel('Im\langleG_2_2(L,j\omega)\rangle')
%title('Nyquist plot for G_2_2(j\omega)')

