function [gttw_hat,gtsw_hat,gstw_hat,gssw_hat] = calc_exchanger_approx_freq_resp(omega,ln,v,k)
%
% This function calculates the spatially distributed frequency responses 
% for the approximate rational transfer function model of the double-pipe heat exchanger 
% for the parallel- (vt>0,vs>0) and counter-flow (vt>0,vs<0) configurations.
% The approximation model consists of N nodes located inside the spatial domain (0,L)
% represented by rational transfer functions evaluated at these nodes.
%
% Function inputs: 
% omega - vector of the angular frequencies (for omega=0 we obtain the steady-state profiles)
% ln - vector of N+2(minimum 3) spatial positions of the MOL approximation model, ln=[0 l1 l2 ... lN L],
% v - vector of the fluid velocities in PDEs, v = [vt vs],
% k - vector of the constant parameters in PDEs, k = [k1 k2 k3 k4].
%
% Function outputs: 
% matrices of the approximate frequency responses evaluated at the nth
% spatial node
% for the parallel flow:
% gttw_hat(ln,i*omega) = Tt_hat(ln,i*omega)/Tt(0,i*omega)
% gtsw_hat(ln,i*omega) = Tt_hat(ln,i*omega)/Ts(0,i*omega)
% gstw_hat(ln,i*omega) = Ts_hat(ln,i*omega)/Tt(0,i*omega)
% gssw_hat(ln,i*omega) = Ts_hat(ln,i*omega)/Ts(0,i*omega)
% and for the counter flow:
% gttw_hat(ln,i*omega) = Tt_hat(ln,i*omega)/Tt(0,i*omega)
% gtsw_hat(ln,i*omega) = Tt_hat(ln,i*omega)/Ts(L,i*omega)
% gstw_hat(ln,i*omega) = Ts_hat(ln,i*omega)/Tt(0,i*omega)
% gssw_hat(ln,i*omega) = Ts_hat(ln,i*omega)/Ts(L,i*omega)%
%
vt = v(1);  vs = abs(v(2));
k1  = k(1); k2  = k(2); k3  = k(3); k4  = k(4);
L = ln(end);  
N = length(ln)-2;

for n=2:N+1                % calculating parameter values for the following sections  
                           % for n=2 we have the first section and for n=N+1 - the last one
    dln = ln(n)-ln(n-1);   % current section length for v(1)>0 and v(2)>0
    if v(2)>0
      dln1 = dln;
    else
      dln1 = ln(n+1)-ln(n);  % current section length for v(2)<0
    end    
    %
    % calculation of the transfer function matrix for the single section
    %
    % coefficients of the denominator polynomial 
    a2 = k1 + k2 + k3 + k4 + vt/dln + vs/dln;
    a1 = k1*k3 + k1*k4 + k2*k4 + vt/dln*(k2 + k3 +k4) + vs/dln1*(k1 + k2 +k3) + vs*vt/dln/dln1;
    a0 = vt/dln*k2*k4 + vs/dln1*k1*k3 + vs*vt*(k2 + k3)/dln/dln1;
    %
    % coefficients of numerator polynomials
    btt_2 = vt/dln;   
    btt_1 = vt*(k2 + k3 + k4)/dln + vs*vt/dln/dln1;  
    btt_0 = vt*k2*k4/dln + vs*vt*(k2 + k3)/dln/dln1;
    bts_0 = vs*k1*k3/dln1; 
    bst_0 = vt*k2*k4/dln;   
    bss_2 = vs/dln1;   
    bss_1 = vs*(k1 + k2 + k3)/dln1 + vs*vt/dln/dln1;  
    bss_0 = vs*k1*k3/dln1 + vs*vt*(k2 + k3)/dln/dln1;
    %
    % coefficient vector of the denominator
    m   = [1 a2 a1 a0];
    % coefficient vector of the numerator
    ltt = [btt_2 btt_1 btt_0];
    lts = [bts_0];
    lst = [bst_0];
    lss = [bss_2 bss_1 bss_0];
    %
    % transfer functions for the single section
    gttn = tf(ltt,m);
    gtsn = tf(lts,m);
    gstn = tf(lst,m);
    gssn = tf(lss,m);
    %
    % transfer function matrix of the single nth section
    Gn = [gttn gtsn; 
          gstn gssn];
    %
    if n==2         % first section 
        Gi = Gn;
    else
        Gi = append(Gi,Gn); % appending of subsequent sections 
    end    
end
%
if v(1)>0 && v(2)>0 % collocated boundary conditions (parallel-flow)
%
% connection matrix for collocated configuration
connections = [ [3:2+2*(N-1)]' [1:2*(N-1)]' ];
% input vector for collocated configuration
inputs = [1 2];
% boundary output vector for collocated configuration
%outputs = [2*N-1 2*N];
% distributed outputs for collocated configuration
outputs = [1:2*N];
%
% section connection for collocated configuration
G_hat = connect(Gi,connections,inputs,outputs);
%
%
elseif v(1)>0 && v(2)<0 % anti-collocated boundary conditions (counter-flow)
%
% connection matrix for anti-collocated configuration
connections = zeros(2*(N-1),2);
connections(1:2:2*N-3,:) = [ [3:2:2*N-1]' [1:2:2*N-3]' ];
connections(2:2:2*N-2,:) = [ [2:2:2*N-2]' [4:2:2*N]' ];
%
% input vector for anti-collocated configuration
inputs = [1 2*N];
%
% boundary output vector for anti-collocated configuration
% outputs = [2*N-1 2];
% distributd outputs for anti-collocated configuration
%outputs = zeros(1,2*N);
%outputs(1,1:2:2*N-1) = [2*N-1:-2:1];
%outputs(1,2:2:2*N) = [2*N:-2:2];
outputs = [1:2*N];
%
% section connection for anti-collocated configuration
G_hat = connect(Gi,connections,inputs,outputs);
%
end
%
% calculation of approximate frequency responses 
F_hat = freqresp(G_hat,omega);
C = permute(F_hat,[1 3 2]);
C = reshape(C,[],size(F_hat,2),1);

gttw_hat = reshape(C(1:2:2*N*length(omega)-1,1),[N length(omega)]).';
gstw_hat = reshape(C(2:2:2*N*length(omega),1),[N length(omega)]).';
gtsw_hat = reshape(C(1:2:2*N*length(omega)-1,2),[N length(omega)]).';
gssw_hat = reshape(C(2:2:2*N*length(omega),2),[N length(omega)]).';

