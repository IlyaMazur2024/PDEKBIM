% ===========================================================================
%
% Skrypt rysuj¹cy charakterystyki czêstotliwoœciowe 2D Bodego i Nyquista 
% uk³adu hiperbolicznego z jedn¹ zmienn¹ przestrzenn¹ w oparciu o model 
% opisany transmitancjami operatorowymi dla uk³adu o parametrach
% roz³o¿onych, warunki pocz¹tkowe zgodne

close all
clear all

init_hyperb                    % inicjalizacja wartoœci modelu 
omega = logspace(-4,1,10000);  % wektor pulsacji dla wymiennika
l = L;                         % charakterystyki dla punktu koñcowego 

% calculation of frequency responses for distributed parameter model
[G11w,G12w,G21w,G22w] = calc_hyperb_distr_freq_tf(omega,l,L,LAMBDA,K);
% calculation of transfer function matrix of the approximation model
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
% calculation of frequency responses for approximation model
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);

% Charakterystyki Bodego

figure(1)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11w)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on 
%title('Amplitude plot for G_1_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G11w)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_1_1(j\omega)')
grid on
hold on 
axis([-4 1 -60 10])
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11Nw)),'b')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G11Nw)),'b')

figure(2)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12w)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on 
%title('Amplitude plot for G_1_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12w)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on 
%title('Phase plot for G_1_2(j\omega)')
axis([-4 1 -60 10])
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12Nw)),'b')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12Nw)),'b')

figure(3)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21w)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on 
%title('Amplitude plot for G_2_1(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21w)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_1(j\omega)')
grid on
hold on 
axis([-4 1 -60 10])
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21Nw)),'b')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21Nw)),'b')

figure(4)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22w)),'r')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on 
%title('Amplitude plot for G_2_2(j\omega)')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22w)),'r')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Phase plot for G_2_2(j\omega)')
grid on
hold on 
axis([-4 1 -60 10])
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22Nw)),'b')
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22Nw)),'b')


% Charakterystyki amplitudowo-fazowe Nyquista
%
figure(5)
plot(real(G11w),imag(G11w),'r');
xlabel('Re\langleG_1_1(L,j\omega)\rangle')
ylabel('Im\langleG_1_1(L,j\omega)\rangle')
%title('Nyquist plot for G_1_1(j\omega)')
grid on
hold on 
plot(real(G11Nw),imag(G11Nw),'b');

figure(6)
plot(real(G12w),imag(G12w),'r');
xlabel('Re\langleG_1_2(L,j\omega)\rangle')
ylabel('Im\langleG_1_2(L,j\omega)\rangle')
%title('Nyquist plot for G_1_2(j\omega)')
grid on
hold on 
plot(real(G12Nw),imag(G12Nw),'b');

figure(7)
plot(real(G21w),imag(G21w),'r');
xlabel('Re\langleG_2_1(L,j\omega)\rangle')
ylabel('Im\langleG_2_1(L,j\omega)\rangle')
%title('Nyquist plot for G_2_1(j\omega)')
grid on
hold on 
plot(real(G21Nw),imag(G21Nw),'b');

figure(8)
plot(real(G22w),imag(G22w),'r');
xlabel('Re\langleG_2_2(L,j\omega)\rangle')
ylabel('Im\langleG_2_2(L,j\omega)\rangle')
%title('Nyquist plot for G_2_2(j\omega)')
grid on
hold on 
plot(real(G22Nw),imag(G22Nw),'b');



