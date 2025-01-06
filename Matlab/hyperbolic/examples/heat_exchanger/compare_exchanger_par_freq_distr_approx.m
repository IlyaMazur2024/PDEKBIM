% ===========================================================================
%
% Comaprison of the frequency responses obtained from the PDE-based 
% irrational transfer function model of the parallel-flow heat exchanger 
% and the approximate rational transfer function models of different orders
% (for different numbers of spatial sections)
%
% Krzysztof Bartecki, 2020

close all
clear variables

omega = logspace(-4,1,10000);   % vector of the angular frequencies
init_exchanger;                 % initialization of the PDE model parameters
l = L;                          % only boundary responses are considered

% calculation of the original frequency responses 
% from the PDE-based irrational transfer function model
[gttw,gtsw,gstw,gssw] = calc_exchanger_par_distr_freq_resp(omega,l,v,k);

% Nyquist plots of the PDE-based frequency responses 

set(0,'DefaultLineLineWidth',1)
set(0,'defaultAxesFontSize',14)

figure(1)
plot(real(gttw),imag(gttw),'k-');
xlabel('Re\langle g_t_t(L,i\omega)\rangle')
ylabel('Im\langle g_t_t(L,i\omega)\rangle')
%title('Nyquist plots for g_t_t(L,i\omega)')
axis([-1 1 -1 0.9])
grid on
hold on

figure(2)
plot(real(gtsw),imag(gtsw),'k-');
%xlabel('Re\langle g_t_s(L,i\omega)\rangle')
xlabel('Re\langle G^{12,u}(i\omega)\rangle')
%ylabel('Im\langle g_t_s(L,i\omega)\rangle')
ylabel('Im\langle G^{12,u}(i\omega)\rangle')
%title('Nyquist plots for g_t_s(L,i\omega)')
axis([-0.07 0.14 -0.12 0.025])
grid on
hold on

figure(3)
plot(real(gstw),imag(gstw),'k-');
xlabel('Re\langle g_s_t(L,i\omega)\rangle')
ylabel('Im\langle g_s_t(L,i\omega)\rangle')
%title('Nyquist plots for g_s_t(L,i\omega)')
axis([-0.25 0.45 -0.35 0.1])
grid on
hold on

figure(4)
plot(real(gssw),imag(gssw),'k-');
xlabel('Re\langle g_s_s(L,i\omega)\rangle')
ylabel('Im\langle g_s_s(L,i\omega)\rangle')
%title('Nyquist plots for g_s_s(L,i\omega)')
axis([-0.7 0.7 -0.7 0.7])
grid on
hold on

% Bode plots of the PDE-based frequency responses 

figure(5)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(gttw)),'k-')
xlabel('lg(\omega)')
ylabel('L_m(L,\omega) [dB]')
%title('Bode magnitude plots for g_t_t(L,i\omega)')
grid on
hold on 
axis([-4 1 -140 10])
subplot(2,1,2)
plot(log10(omega),unwrap(angle(gttw)),'k-')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
title('Bode phase plots for g_t_t(L,i\omega)')
grid on
hold on 
axis([-4 1 -50 10])

figure(6)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(gtsw)),'k-')
xlabel('lg(\omega)')
ylabel('L_m(L,\omega) [dB]')
%title('Bode magnitude plots for g_t_s(L,i\omega)')
grid on
hold on 
axis([-4 1 -140 10])
subplot(2,1,2)
plot(log10(omega),unwrap(angle(gtsw)),'k-')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Bode phase plots for g_t_s(L,i\omega)')
grid on
hold on 
axis([-4 1 -50 10])

figure(7)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(gstw)),'k-')
xlabel('lg(\omega)')
ylabel('L_m(L,\omega) [dB]')
%title('Bode magnitude plots for g_s_t(L,i\omega)')
grid on
hold on 
axis([-4 1 -140 10])
subplot(2,1,2)
plot(log10(omega),unwrap(angle(gstw)),'k-')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Bode phase plots for g_s_t(L,i\omega)')
grid on
hold on 
axis([-4 1 -50 10])

figure(8)
subplot(2,1,1)
plot(log10(omega),20*log10(abs(gssw)),'k-')
xlabel('lg(\omega)')
ylabel('L_m(L,\omega) [dB]')
%title('Bode magnitude plots for g_s_s(L,i\omega)')
grid on
hold on 
axis([-4 1 -140 10])
subplot(2,1,2)
plot(log10(omega),unwrap(angle(gssw)),'k-')
xlabel('lg(\omega)')
ylabel('\phi(L,\omega) [rad]')
%title('Bode phase plots for g_s_s(L,i\omega)')
grid on
hold on 
axis([-4 1 -50 10])

% calculation of the approximate frequency responses 
% from rational transfer function models of different orders

% Four different numbers of uniform sections 
for N=[1 9 99 999]
%     dl = L/N;
%     ln = [dl:dl:L]; % vector of N (uniform) spatial positions of the MOL approximation model
%     [gttw_hat,gtsw_hat,gstw_hat,gssw_hat] = calc_exchanger_par_approx_freq_resp(omega,ln,v,k);
%     % only boundary responses are considered
%     gttw_hat = gttw_hat(:,end);
%     gtsw_hat = gtsw_hat(:,end);
%     gstw_hat = gstw_hat(:,end);
%     gssw_hat = gssw_hat(:,end);
     
    % plot markers and colors for different N
   switch N    
        case 1
          load G1w_exchanger_par
           s = 'b:';
        case 9
          load G10w_exchanger_par
           s = 'r-.';
        case 99
          load G100w_exchanger_par
           s = 'g--';
        case 999
          load G1000w_exchanger_par
           s = 'm:';
    end

    %  Nyquist plots of the approximate frequency responses 
    figure(1)
    plot(real(gttw_hat),imag(gttw_hat),s);
    figure(2)
    plot(real(gtsw_hat),imag(gtsw_hat),s);
    figure(3)
    plot(real(gstw_hat),imag(gstw_hat),s);
    figure(4)
    plot(real(gssw_hat),imag(gssw_hat),s);

    %  Bode plots of the approximate frequency responses 
    figure(5)
    subplot(2,1,1)
    plot(log10(omega),20*log10(abs(gttw_hat)),s)
    subplot(2,1,2)
    plot(log10(omega),unwrap(angle(gttw_hat)),s)
    figure(6)
    subplot(2,1,1)
    plot(log10(omega),20*log10(abs(gtsw_hat)),s)
    subplot(2,1,2)
    plot(log10(omega),unwrap(angle(gtsw_hat)),s)
    figure(7)
    subplot(2,1,1)
    plot(log10(omega),20*log10(abs(gstw_hat)),s)
    subplot(2,1,2)
    plot(log10(omega),unwrap(angle(gstw_hat)),s)
    figure(8)
    subplot(2,1,1)
    plot(log10(omega),20*log10(abs(gssw_hat)),s)
    subplot(2,1,2)
    plot(log10(omega),unwrap(angle(gssw_hat)),s)

end

% for i=1:4
%     figure(i)
%     legend('PDE','N=1','N=10','N=100','Location','southwest','Orientation','horizontal');    
% end
% for i=5:8
%     figure(i)
%     subplot(2,1,1)
%     legend('PDE','N=1','N=10','N=100','Location','southwest','Orientation','horizontal');    
%     subplot(2,1,2)
%     legend('PDE','N=1','N=10','N=100','Location','southwest','Orientation','horizontal');    
% end
for i=1:8
    figure(i)
    legend('PDE','N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal');    
    %legend('PDE','N=1','N=10','N=100','Location','northoutside','Orientation','horizontal');    
end
