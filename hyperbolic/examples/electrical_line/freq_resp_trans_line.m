function [G11c,G12c,G21c,G22c,G11i,G12i,G21i,G22i] = freq_resp_trans_line(omega,l,R,L,G,C)



Lp = l(end)*ones(1,length(l));

gamma = sqrt((R+j*omega*L).*(G+j*omega*C));
Z = sqrt((R+j*omega*L)./(G+j*omega*C));   
Y = 1./Z;   

G11c = cosh(gamma.'*l);
G12c = -repmat(Y.',1,length(l)).*sinh(gamma.'*l);
G21c = -repmat(Z.',1,length(l)).*sinh(gamma.'*l);
G22c = cosh(gamma.'*Lp);

G11i = (exp(2*gamma.'*Lp)+exp(2*gamma.'*l))./(exp(gamma.'*l).*(exp(2*gamma.'*Lp)+1));
G12i = -repmat(Z.',1,length(l)).*(sinh(gamma.'*l)./cosh(gamma.'*Lp)) ;
G21i =  repmat(Y.',1,length(l)).*(exp(2*gamma.'*Lp)-exp(2*gamma.'*l))./(exp(gamma.'*l).*(exp(2*gamma.'*Lp)+1));
G22i = cosh(gamma.'*l)./cosh(gamma.'*Lp);
