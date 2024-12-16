function [gttw_hat,gtsw_hat,gstw_hat,gssw_hat] = calc_exchanger_par_approx_freq_resp(omega,ln,v,k)
%
% This function calculates spatially distributed frequency responses 
% for the approximate rational transfer function model of the 
% double-pipe parallel-flow heat exchanger.  
% The approximation model consists of N spatial (possibly non-uniform) sections. 
%
% Function inputs: 
% omega - vector of the angular frequencies, 
% ln - vector of N spatial positions of the approximation model, ln=[l1 l2 ... lN],
% v - vector of the fluid velocities in PDEs, v = [vt vs],
% k - vector of the constant parameters in PDEs, k = [k1 k2 k3 k4].
%
% Function outputs: 
% matrices of the approximate frequency responses evaluated 
% for each nth section for the parallel-flow configuration:
% gttw_hat(ln,i*omega) = Tt_hat(ln,i*omega)/Tt(0,i*omega)
% gtsw_hat(ln,i*omega) = Tt_hat(ln,i*omega)/Ts(0,i*omega)
% gstw_hat(ln,i*omega) = Ts_hat(ln,i*omega)/Tt(0,i*omega)
% gssw_hat(ln,i*omega) = Ts_hat(ln,i*omega)/Ts(0,i*omega)
%
% Krzysztof Bartecki, 2020

vt = v(1);  vs = v(2);
k1  = k(1); k2  = k(2); k3  = k(3); k4  = k(4);
L = ln(end);  
N = length(ln);

for n=1:N
    if n==1 
        dl=ln(1);              % section length for n=1
    else    
        dl = ln(n)-ln(n-1);    % section lengths for n=2,3,...,N 
    end    
    
    % calculation of the transfer function matrix for the single section
    
    % coefficients of the denominator polynomial 
    a2 = k1 + k2 + k3 + k4 + (vt + vs)/dl;
    a1 = k1*k3 + k1*k4 + k2*k4 + (vt*(k2 + k3 +k4) + vs*(k1 + k2 +k3))/dl + vs*vt/(dl^2);
    a0 = (vt*k2*k4 + vs*k1*k3)/dl + vs*vt*(k2 + k3)/(dl^2);
    
    % coefficients of numerator polynomials
    btt_2 = vt/dl;   
    btt_1 = vt*(k2 + k3 + k4)/dl + vs*vt/(dl^2);  
    btt_0 = vt*k2*k4/dl + vs*vt*(k2 + k3)/(dl^2);
    bts_0 = vs*k1*k3/dl; 
    bst_0 = vt*k2*k4/dl;   
    bss_2 = vs/dl;   
    bss_1 = vs*(k1 + k2 + k3)/dl + vs*vt/(dl^2);  
    bss_0 = vs*k1*k3/dl + vs*vt*(k2 + k3)/(dl^2);
    
    % coefficient vector of the denominator
    m   = [1 a2 a1 a0];
    % coefficient vector of the numerator
    ltt = [btt_2 btt_1 btt_0];
    lts = [bts_0];
    lst = [bst_0];
    lss = [bss_2 bss_1 bss_0];
    
    % transfer functions for the single section
    gttn = tf(ltt,m);
    gtsn = tf(lts,m);
    gstn = tf(lst,m);
    gssn = tf(lss,m);
    
    % transfer function matrix of the single nth section
    Gn = [gttn gtsn; 
          gstn gssn];
    %
    if n==1 
        Gi = Gn;
    else
        Gi = append(Gi,Gn); % appending of subsequent sections 
    end    
end

if vt>0 && vs>0      % parallel-flow configuration

% connection matrix for collocated parallel-flow 
connections = [ [3:2+2*(N-1)]' [1:2*(N-1)]' ];

% input vector for parallel-flow 
inputs = [1 2];

% boundary output vector for parallel-flow 
% outputs = [2*N-1 2*N];
% distributed outputs for parallel-flow 
outputs = [1:2*N];

% section connection for parallel-flow 
G_hat = connect(Gi,connections,inputs,outputs);

% calculation of approximate frequency responses 
F_hat = freqresp(G_hat,omega);
C = permute(F_hat,[1 3 2]);
C = reshape(C,[],size(F_hat,2),1);

gttw_hat = reshape(C(1:2:2*N*length(omega)-1,1),[N length(omega)]).';
gstw_hat = reshape(C(2:2:2*N*length(omega),1),[N length(omega)]).';
gtsw_hat = reshape(C(1:2:2*N*length(omega)-1,2),[N length(omega)]).';
gssw_hat = reshape(C(2:2:2*N*length(omega),2),[N length(omega)]).';

else
    disp('This function works only for vt>0 and vs>0.')
end

