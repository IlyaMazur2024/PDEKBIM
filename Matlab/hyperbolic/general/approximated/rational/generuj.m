clear all
init_hyperb

omega = logspace(-4,2,10000);
T = 100;
delta_t=0.01;
t = [0:delta_t:T];  
l = L;              

lambda2 = 0.2;
LAMBDA = diag([lambda1, lambda2]);

N = 1;
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);
[g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);
save G1_col.mat 

N = 10;
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);
[g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);

save G10_col.mat 

N = 100;
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);
[g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);

save G100_col.mat 

N = 1000;
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);
[g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);

save G1000_col.mat 


lambda2 = -0.2;
LAMBDA = diag([lambda1, lambda2]);

N = 1;
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);
[g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);
save G1_ant.mat 

N = 10;
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);
[g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);

save G10_ant.mat 

N = 100;
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);
[g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);

save G100_ant.mat 

N = 1000;
[GN] = calc_hyperb_approx_tf(LAMBDA,K,L,N);
[G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyperb_approx_freq_tf(GN,omega);
[g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t);

save G1000_ant.mat 

