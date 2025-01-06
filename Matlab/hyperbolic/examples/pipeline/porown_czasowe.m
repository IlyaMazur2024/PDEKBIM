% ===========================================================================
%
% Skrypt porównuj¹cy odpowiedzi czasowe ruroci¹gu sprawnego i uszkodzonego
% dla transmitancji G(s,L) = P(s,L)/Q(s,L) = - Z*tanh(gamma*L)
% 
close all
clear all
clc

Lp=1000;
omega = logspace(-4,2,10000000);
t = [0:0.1:30];
P = 5E5;           % ciœnienie na wlocie
q = 27.5;          % przep³yw w kg/s

G=trans_widm(omega,Lp,q);  % transmitacja widmowa ruroci¹gu bez uszkodzeñ

load odp_skok

% % odpowiedŸ skokowa ruroci¹gu bez uszkodzeñ
% hG1 = zeros(4,length(t));
% % obliczanie wartoœci wyra¿enia podca³kowego dla kolejnych chwil czasowych
% for i = 1:length(t)
%    wart_hG = P_sin_omega(G,omega,t(i));
%    % ca³kowanie po omega  
%    hG1(1,i) = trapz(omega,wart_hG);
% end
% 
% 
% % odpowiedzi skokowe ruroci¹gu z uszkodzeniami dla ró¿nych œrednic wycieku
% l1 = 750;                % pozycja wycieku
% Dr = [0.01 0.02 0.1];    % œrednice wycieku
% 
% for i=1:3
%     i
%    Dl=Dr(i);
%    % oblicza spadki ciœnieñ, strumienie oraz rezystancje hydrauliczne w stanie ustalonym
%    [dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,Lp,l1,Dl);
%    % transmitacja widmowa ruroci¹gu z wyciekiem
%    Gl=trans_widm_uszk(omega,l1,Dl,q1,q2,ql);
%       
%    for j = 1:length(t)
%      j
%      wart_hG = P_sin_omega(Gl,omega,t(j));
%      % ca³kowanie po omega  
%      hG1(i+1,j) = trapz(omega,wart_hG);
%    end
% end
% 
% hG1 = 2/pi*hG1;


% wykres
%col = ['k' 'r' 'c' 'b'];
col = {'k-' 'k:' 'k-.' 'k--'};

figure(1)
for i=1:4
   plot(t,-hG1(i,:),col{i})
   hold on
end
xlabel('t(s)')
ylabel('p(L_p,t) [Pa]')


% odpowiedzi skokowe ruroci¹gu z uszkodzeniami dla ró¿nych lokalizacji wycieku

% hG2 = zeros(4,length(t));


% Dl = 0.02;   % Dl=D/10
% lr = [250 500 750];
% 
% for i=1:3
%     i
%    l1=lr(i);  
%   
%    % oblicza spadki ciœnieñ, strumienie oraz rezystancje hydrauliczne w stanie ustalonym
%    [dp1,dp2,dpl,dpk,q1,q2,ql,R1,R2,Rl,Rk,Re1,Re2,Rel,Rek] = oblicz_ustalony(P,Lp,l1,Dl);
%    % transmitacja widmowa ruroci¹gu z wyciekiem
%    Gl=trans_widm_uszk(omega,l1,Dl,q1,q2,ql);
%    
%    for j = 1:length(t)
%      j
%      wart_hG = P_sin_omega(Gl,omega,t(j));
%      % ca³kowanie po omega  
%      hG2(i+1,j) = trapz(omega,wart_hG);
%    end
% end
% 
% hG2 = 2/pi*hG2;
% hG2(1,:) = hG1(1,:);

% wykres
%col = ['k' 'r' 'c' 'b'];
col = {'k-' 'k:' 'k-.' 'k--'};

figure(2)
for i=1:4
   plot(t,-hG2(i,:),col{i})
   hold on
end
xlabel('t(s)')
ylabel('p(L_p,t) [Pa]')





