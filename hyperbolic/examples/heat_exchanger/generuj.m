clear all
init_exchanger
omega = 0;
%logspace(-4,1,10000);
        
N = 1;
ln = [0:L/(N+1):L];
[gttw_hat,gtsw_hat,gstw_hat,gssw_hat] = calc_exchanger_approx_freq_resp(omega,ln,v,k);
save G1w_exchanger_par.mat 

N = 9;
ln = [0:L/(N+1):L];
[gttw_hat,gtsw_hat,gstw_hat,gssw_hat] = calc_exchanger_approx_freq_resp(omega,ln,v,k);
save G10w_exchanger_par.mat 

N = 99;
ln = [0:L/(N+1):L];
[gttw_hat,gtsw_hat,gstw_hat,gssw_hat] = calc_exchanger_approx_freq_resp(omega,ln,v,k);
save G100w_exchanger_par.mat 

 N = 999;
 ln = [0:L/(N+1):L];
 [gttw_hat,gtsw_hat,gstw_hat,gssw_hat] = calc_exchanger_approx_freq_resp(omega,ln,v,k);
 save G1000w_exchanger_par.mat 