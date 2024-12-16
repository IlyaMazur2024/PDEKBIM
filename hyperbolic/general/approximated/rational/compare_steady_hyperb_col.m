% ===========================================================================
%
% Skrypt porównuj¹cy profile zmiennych stanu uk³adu hiperbolicznego 2x2 
% dla stanu ustalonego, dla uk³adu o parametrach roz³o¿onych oraz modeli 
% aproiksymacyjnych , warunki pocz¹tkowe zgodne

close all
clear all
set(0,'DefaultLineLineWidth',1.5)
set(0,'defaultAxesFontSize',14)

init_hyperb                    % inicjalizacja wartoœci modelu 
omega = 0;                     % stan ustalony
l = [0:0.1:L];                 % charakterystyki dla punktu koñcowego 

x10 = 100;   % warunki brzegowe (temperatury na wlocie)
x20 = 50;     

% approximation models for N=1, 10, 100, 1000
H = gobjects(1,5);
i = 0;

for N = [1 10 100 1000]
  switch N    
    case 1
        s = 'b:';
    case 10
        s = 'r-.';
    case 100
        s = 'g--';
    case 1000
        s = 'm:';
  end

dl = L/N;
% calculation of transfer function matrix of the approximation model
[GN] = calc_hyperb_approx_ss(LAMBDA,K,L/N,1);
%GN=minreal(GN,0.01);
% calculation of frequency responses for approximation model
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);

x1Ns = x10;
x2Ns = x20;
dl = L/N;
l = 0;

for n=1:N
    x1Ns_a = G11Nw*x1Ns + G12Nw*x2Ns; 
    x2Ns_a = G21Nw*x1Ns + G22Nw*x2Ns; 
    
    x1Ns = x1Ns_a; 
    x2Ns = x2Ns_a; 
    
    figure(1)
    g=plot([l l+dl],[x1Ns x1Ns],s,'LineWidth',1.5);
    if N==100
        g.Color = [0 0.5 0];
    end
    hold on
    h=plot([l l+dl],[x2Ns x2Ns],s,'LineWidth',1.5);    
    if N==100
        h.Color = [0 0.5 0];
    end
    l=l+dl;
end
i = i+1;    
H(i)= h;
end

% calculation of frequency responses for distributed parameter model
l = [0:0.1:L];
[G11w,G12w,G21w,G22w] = calc_hyperb_distr_freq_tf(omega,l,L,LAMBDA,K);
% steady state profiles for distributed parameter model
x1s = G11w*x10 + G12w*x20; 
x2s = G21w*x10 + G22w*x20; 

figure(1)
plot(l,x1s,'k-')
hold on
h=plot(l,x2s,'k-');
H(5)=h;
xlabel('l [m]','FontSize',14)
ylabel('\vartheta_1(l),  \vartheta_2(l)','FontSize',14)
grid on
text(0.5,52,'\vartheta_2(l)','FontSize',14)
text(0.5,102,'\vartheta_1(l)','FontSize',14)
axis([0 5 48 105])
legend(H,'N=1','N=10','N=100','N=1000','PDE','Location','northoutside','Orientation','horizontal')


