% ===========================================================================
%
% Skrypt porównuj¹cy profile zmiennych stanu uk³adu hiperbolicznego 2x2 
% dla stanu ustalonego, dla uk³adu o parametrach roz³o¿onych oraz modeli 
% aproksymacyjnych, warunki brzegowe przeciwne

close all
clear all

init_hyperb                    % inicjalizacja wartoœci modelu 
omega = 0;                     % stan ustalony
l = [0:0.1:L];                 % charakterystyki dla punktu koñcowego 

x10 = 100;   % warunki brzegowe (temperatury na wlocie)
x2L = 50;     

% approximation models for N=1, 10, 100, 1000
H = gobjects(1,5);
i = 0;

for N = [10 100 1000]

  switch N    
    case 1
        s = 'b:';
    case 10
        s = 'r-.';
    case 100
        s = 'g--';
    case 1000
        s = 'm-';
  end


% calculation of transfer function matrix of the approximation model
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
% calculation of frequency responses for approximation model
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);

u=[x10; x2L];
StaticGain=dcgain(GN);
y=StaticGain*u;

%x1Ns=y(1:2:2*N);
%x2Ns=y(2:2:2*N);

x1Ns=y(1);
x2Ns=y(2);

l = 0;
dl = L/N;

figure(1)
for n=1:N
    
  
    g=plot([l l+dl],[x1Ns(n) x1Ns(n)],s,'LineWidth',1.5);
    if N==100
        g.Color = [0 0.5 0];
    end
    hold on
    h=plot([l l+dl],[x2Ns(n) x2Ns(n)],s,'LineWidth',1.5);    
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
x1s = G11w*x10 + G12w*x2L; 
x2s = G21w*x10 + G22w*x2L; 

figure(1)
plot(l,x1s,'k-')
hold on
h=plot(l,x2s,'k-');
H(5)=h;
xlabel('l [m]','FontSize',14)
ylabel('\vartheta_1(l),  \vartheta_2(l)','FontSize',14)
grid on
text(0.5,77,'\vartheta_2(l)','FontSize',12)
text(0.5,102,'\vartheta_1(l)','FontSize',12)
axis([0 5 48 105])
legend(H,'N=1','N=10','N=100','N=1000','DPS','Location','northoutside','Orientation','horizontal')


