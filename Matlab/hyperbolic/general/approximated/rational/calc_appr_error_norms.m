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

load G1_col
s = 'k:';

e11_1 = G11w - G11Nw;
figure(1)        
plot(log10(omega),abs(e11_1),s)
grid on
hold on
e_2_11_1 = sqrt(1/(2*pi)*sum(abs(e11_1).^2))
[e_inf_11_1 i_inf_11_1] = max(abs(e11_1))        
        
e12_1 = G12w - G12Nw;
figure(2)        
plot(log10(omega),abs(e12_1),s)
grid on
hold on
e_2_12_1 = sqrt(1/(2*pi)*sum(abs(e12_1).^2))
[e_inf_12_1 i_inf_12_1] = max(abs(e12_1))
        
e21_1 = G21w - G21Nw;
figure(3)        
plot(log10(omega),abs(e21_1),s)
grid on
hold on
e_2_21_1 = sqrt(1/(2*pi)*sum(abs(e21_1).^2))
[e_inf_21_1 i_inf_21_1] = max(abs(e21_1))
       
e22_1 = G22w - G22Nw;
figure(4)        
plot(log10(omega),abs(e22_1),s)
grid on
hold on
e_2_22_1 = sqrt(1/(2*pi)*sum(abs(e22_1).^2))
[e_inf_22_1 i_inf_22_1] = max(abs(e22_1))        


load G10_col
s = 'b-.';

e11_10 = G11w - G11Nw;
figure(1)        
plot(log10(omega),abs(e11_10),s)
e_2_11_10 = sqrt(1/(2*pi)*sum(abs(e11_10).^2))
[e_inf_11_10 i_inf_11_10] = max(abs(e11_10))   

e12_10 = G12w - G12Nw;
figure(2)        
plot(log10(omega),abs(e12_10),s)
e_2_12_10 = sqrt(1/(2*pi)*sum(abs(e12_10).^2))
[e_inf_12_10 i_inf_12_10] = max(abs(e12_10))

e21_10 = G21w - G21Nw;
figure(3)        
plot(log10(omega),abs(e21_10),s)
e_2_21_10 = sqrt(1/(2*pi)*sum(abs(e21_10).^2))
[e_inf_21_10 i_inf_21_10] = max(abs(e21_10))

e22_10 = G22w - G22Nw;
figure(4)        
plot(log10(omega),abs(e22_10),s)
e_2_22_10 = sqrt(1/(2*pi)*sum(abs(e22_10).^2))
[e_inf_22_10 i_inf_22_10] = max(abs(e22_10))    

load G100_col
s = 'r--';

e11_100 = G11w - G11Nw;
figure(1)        
plot(log10(omega),abs(e11_100),s)
e_2_11_100 = sqrt(1/(2*pi)*sum(abs(e11_100).^2))
[e_inf_11_100 i_inf_11_100] = max(abs(e11_100))   

e12_100 = G12w - G12Nw;
figure(2)        
plot(log10(omega),abs(e12_100),s)
e_2_12_100 = sqrt(1/(2*pi)*sum(abs(e12_100).^2))
[e_inf_12_100 i_inf_12_100] = max(abs(e12_100))

e21_100 = G21w - G21Nw;
figure(3)        
plot(log10(omega),abs(e21_100),s)
e_2_21_100 = sqrt(1/(2*pi)*sum(abs(e21_100).^2))
[e_inf_21_100 i_inf_21_100] = max(abs(e21_100))

e22_100 = G22w - G22Nw;
figure(4)        
plot(log10(omega),abs(e22_100),s)
e_2_22_100 = sqrt(1/(2*pi)*sum(abs(e22_100).^2))
[e_inf_22_100 i_inf_22_100] = max(abs(e22_100)) 

load G1000_col
s = 'm:';

e11_1000 = G11w - G11Nw;
figure(1)        
plot(log10(omega),abs(e11_1000),s)
axis([-4 2 -0.1 1.1])
xlabel('lg(\omega)')
ylabel('| e_1_1(i\omega) |')
legend('N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')
e_2_11_1000 = sqrt(1/(2*pi)*sum(abs(e11_1000).^2))
[e_inf_11_1000 i_inf_11_1000] = max(abs(e11_1000))   

e12_1000 = G12w - G12Nw;
figure(2)        
plot(log10(omega),abs(e12_1000),s)
axis([-4 2 -0.005 0.06])
xlabel('lg(\omega)')
ylabel('| e_1_2(i\omega) |')
legend('N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')
e_2_12_1000 = sqrt(1/(2*pi)*sum(abs(e12_1000).^2))
[e_inf_12_1000 i_inf_12_1000] = max(abs(e12_1000))

e21_1000 = G21w - G21Nw;
figure(3)        
plot(log10(omega),abs(e21_1000),s)
xlabel('lg(\omega)')
ylabel('| e_2_1(i\omega) |')
legend('N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')
e_2_21_1000 = sqrt(1/(2*pi)*sum(abs(e21_1000).^2))
[e_inf_21_1000 i_inf_21_1000] = max(abs(e21_1000))

e22_1000 = G22w - G22Nw;
figure(4)        
plot(log10(omega),abs(e22_1000),s)
xlabel('lg(\omega)')
ylabel('| e_2_2(i\omega) |')
legend('N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')
e_2_22_1000 = sqrt(1/(2*pi)*sum(abs(e22_1000).^2))
[e_inf_22_1000 i_inf_22_1000] = max(abs(e22_1000)) 
      


