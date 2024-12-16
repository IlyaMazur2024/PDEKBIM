% ===========================================================================
%
% Skrypt rysuj¹cy charakterystyki czêstotliwoœciowe 2D Bodego i Nyquista 
% uk³adu hiperbolicznego z jedn¹ zmienn¹ przestrzenn¹ w oparciu o model 
% opisany transmitancjami operatorowymi dla uk³adu o parametrach
% roz³o¿onych, warunki pocz¹tkowe zgodne

close all
clear all
set(0,'DefaultLineLineWidth',1.5)
set(0,'defaultAxesFontSize',14)

init_hyperb                    % inicjalizacja wartoœci modelu 
omega = logspace(-4,2,10000);  % wektor pulsacji dla wymiennika
l = L;                         % charakterystyki dla punktu koñcowego 

% calculation of frequency responses for distributed parameter model
[G11w,G12w,G21w,G22w] = calc_hyperb_distr_freq_tf(omega,l,L,LAMBDA,K);

% frequency responses for the approximated models, N=1, N=10, N=100, N=1000.
for i=1:4
    
switch i    
    case 1
        load G1_ant
        s = 'b:';
    case 2
        load G10_ant
        s = 'r-.';
    case 3
        load G100_ant
        s = 'g--';
    case 4
        load G1000_ant
        s = 'm:';
end

% Bode plots
figure(1)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11Nw)),s)
hold on 
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G11Nw)),s)
hold on 

figure(2)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12Nw)),s)
hold on 
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12Nw)),s)
hold on 

figure(3)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21Nw)),s)
hold on 
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21Nw)),s)
hold on 

figure(4)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22Nw)),s)
hold on 
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22Nw)),s)
hold on 

% Bode error plots

figure(5)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11w-G11Nw)),s)
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G11w-G11Nw)),s)
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on
legend('N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')

figure(6)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12w-G12Nw)),s)
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12w-G12Nw)),s)
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on
%axis([-4 1 -180 10])
legend('N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')

figure(7)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21w-G21Nw)),s)
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21w-G21Nw)),s)
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on
%axis([-4 1 -180 10])
legend('N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')

figure(8)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22w-G22Nw)),s)
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22w-G22Nw)),s)
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on
axis([-4 1 -290 10])
legend('N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')

% Nyquist plots

figure(9)
plot(real(G11Nw),imag(G11Nw),s);
hold on 

figure(10)
plot(real(G12Nw),imag(G12Nw),s);
hold on 

figure(11)
plot(real(G21Nw),imag(G21Nw),s);
hold on 

figure(12)
plot(real(G22Nw),imag(G22Nw),s);
hold on 

end


% frequency responses for distributed parameter model

% Bode plots

figure(1)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G11w)),'k-')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on 
axis([-4 1 -140 10])
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G11w)),'k-')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on 
axis([-4 1 -50 10])
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')

figure(2)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G12w)),'k-')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on 
axis([-4 1 -100 -10])
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G12w)),'k-')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on 
axis([-4 1 -3.5 0.5])
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')

figure(3)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G21w)),'k-')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on 
axis([-4 1 -80 10])
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G21w)),'k-')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on 
axis([-4 1 -3.5 0.5])
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')

figure(4)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(G22w)),'k-')
xlabel('lg(\omega)')
ylabel('A_v(L,\omega) [dB]')
grid on
hold on 
axis([-4 1 -350 10])
subplot(2,1,2)
plot(log10(omega),unwrap(angle(G22w)),'k-')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
grid on
hold on 
axis([-4 1 -25 10])
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')

% Nyquist plots

figure(9)
plot(real(G11w),imag(G11w),'k-');
xlabel('Re\langleG_1_1(L,i\omega)\rangle','FontSize',14)
ylabel('Im\langleG_1_1(L,i\omega)\rangle','FontSize',14)
grid on
hold on 
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')

figure(10)
plot(real(G12w),imag(G12w),'k-');
xlabel('Re\langleG_1_2(L,i\omega)\rangle','FontSize',14)
ylabel('Im\langleG_1_2(L,i\omega)\rangle','FontSize',14)
grid on
hold on 
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')

figure(11)
plot(real(G21w),imag(G21w),'k-');
xlabel('Re\langleG_2_1(L,i\omega)\rangle','FontSize',14)
ylabel('Im\langleG_2_1(L,i\omega)\rangle','FontSize',14)
grid on
hold on 
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')

figure(12)
plot(real(G22w),imag(G22w),'k-');
xlabel('Re\langleG_2_2(L,i\omega)\rangle','FontSize',14)
ylabel('Im\langleG_2_2(L,i\omega)\rangle','FontSize',14)
grid on
hold on 
legend('N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')


