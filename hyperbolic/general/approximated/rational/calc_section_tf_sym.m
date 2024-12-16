
% skrypt oblicza symboliczne wyra¿enia na transmitancje operatorowe
% pojedynczej sekcji modelu aproksymacyjnego 

syms x1n x1n_1 x2n x2n_1 k11 k12 k21 k22 real
syms dl lam1 lam2 positive
syms a1 a0
syms b11_1 b11_0 b12_0 b21_0 b22_1 b22_0
syms m l11 l12 l21 l22
syms s G11n G12n G21n G22n Gn

a1 = lam1/dl + lam2/dl - k11 - k22;
a0 = (k11 - lam1/dl)*(k22 - lam2/dl)-k21*k12;

b11_1 = lam1/dl;      b11_0 = lam1/dl*(lam2/dl-k22);
b12_0 = k12*lam2/dl;  b21_0 = k21*lam1/dl; 
b22_1 = lam2/dl;      b22_0 = lam2/dl*(lam1/dl-k11);

m   = poly2sym([1 a1 a0],s);
l11 = poly2sym([b11_1 b11_0],s);
l12 = poly2sym([b12_0],s);
l21 = poly2sym([b21_0],s);
l22 = poly2sym([b22_1 b22_0],s);

G11n = l11/m;
G12n = l12/m;
G21n = l21/m;
G22n = l22/m;

Gn = [G11n G12n; G21n G22n];

G10 = Gn^10;

