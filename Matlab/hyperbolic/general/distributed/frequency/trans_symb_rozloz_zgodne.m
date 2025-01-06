
syms k11 k12 k21 k22 lambda1 lambda2 l k real
syms s complex

s1 = 1/2*(k11+k22+sqrt((k11-k22)^2+4*k12*k21));
s2 = 1/2*(k11+k22-sqrt((k11-k22)^2+4*k12*k21));

alpha = -1/(2*lambda1)*(s-k11)-1/(2*lambda2)*(s-k22);
beta = 1/2*sqrt( (1/lambda1*(s-k11)+ 1/lambda2*(s-k22))^2 - 4/(lambda1*lambda2)*(s-s1)*(s-s2) );

q1 = alpha + beta;
q2 = alpha - beta;

A1 = (q1+1/lambda2*(s-k22))/(q1-q2);
B1 = (q2+1/lambda2*(s-k22))/(q1-q2);
A2 = (q1+1/lambda1*(s-k11))/(q1-q2);
B2 = (q2+1/lambda1*(s-k11))/(q1-q2);

G11 = A1* exp(q1*l) - B1 * exp(q2*l);
G12 = k12/lambda1/(q1-q2)*(exp(q1*l)-exp(q2*l));
G21 = k21/lambda2/(q1-q2)*(exp(q1*l)-exp(q2*l));
G22 = A2*exp(q1*l) - B2*exp(q2*l);


G110000 = subs(G11,{k12 k21},{0 0})
G120000 = subs(G12,{k12 k21},{0 0})
G210000 = subs(G21,{k12 k21},{0 0})
G220000 = subs(G22,{k12 k21},{0 0})

% bieguny transmitancji G11, G12, G21 i G22
% oraz zera transmitancji G11 i G22
root1 = (k22*lambda1-k11*lambda2+2*sqrt(-k12*k21*lambda1*lambda2))/(lambda1-lambda2);
root2 = (k22*lambda1-k11*lambda2-2*sqrt(-k12*k21*lambda1*lambda2))/(lambda1-lambda2);

root1n = subs(root1,{lambda1 lambda2 k11 k12 k21 k22},{1 0.2 -0.0638 0.0638 0.0359 -0.0359})
root2n = subs(root2,{lambda1 lambda2 k11 k12 k21 k22},{1 0.2 -0.0638 0.0638 0.0359 -0.0359})

% Zera transmitancji G12 oraz G21
zero1 = 1/(l*(lambda1-lambda2))*((k22*lambda1-k11*lambda2)*l+2*sqrt(-lambda1^2*lambda2^2*k^2*pi^2-lambda1*lambda2*k12*k21*l^2));
zero2 = 1/(l*(lambda1-lambda2))*((k22*lambda1-k11*lambda2)*l-2*sqrt(-lambda1^2*lambda2^2*k^2*pi^2-lambda1*lambda2*k12*k21*l^2));

zero1n = subs(zero1,{lambda1 lambda2 k11 k12 k21 k22 l k},{1 0.2 -0.0638 0.0638 0.0359 -0.0359 5 2})
zero2n = subs(zero2,{lambda1 lambda2 k11 k12 k21 k22 l k},{1 0.2 -0.0638 0.0638 0.0359 -0.0359 5 2})