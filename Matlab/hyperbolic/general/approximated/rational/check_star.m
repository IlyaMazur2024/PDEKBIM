L=0.05; N=1;
GN = calc_hyperb_approx_ss(LAMBDA,K,L,N);
g11 = GN(1,1);
g12 = GN(1,2);
g21 = GN(2,1);
g22 = GN(2,2);

[g11c,g12c,g21c,g22c] = redhefferstar(g11,g12,g21,g22,g11,g12,g21,g22)
Gc = [g11c g12c;g21c g22c];

L=0.1; N=2;
GNcc = calc_hyperb_approx_tf(LAMBDA,K,L,N);
g11cc = GNcc(1,1);
g12cc = GNcc(1,2);
g21cc = GNcc(2,1);
g22cc = GNcc(2,2);

figure(1)
step(g11c)
hold on
step(g22c)
step(g11cc)
step(g22cc)

figure(2)
step(g12c)
hold on
step(g21c)
step(g12cc)
step(g21cc)



