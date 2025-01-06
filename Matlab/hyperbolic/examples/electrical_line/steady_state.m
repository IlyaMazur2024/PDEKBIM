% ===========================================================================
%
% Skrypt porównuj¹cy rozk³ad zmiennych procesowych w stanie ustalonym
% uk³adu hiperbolicznego z jedn¹ zmienn¹ przestrzenn¹ w oparciu o model 
% opisany transmitancjami operatorowymi dla uk³adu o parametrach
% roz³o¿onych, warunki brzegowe zgodne i przeciwstawne
%
clc
close all
clear all

trans_line_init        % inicjalizacja wartoœci modelu linii

% wyznaczenie odpowiedzi uk³adu oryginalnego (sprzê¿onego)

% w oparciu o oryginalne równania 
l = [0:10:Lp];                        
ql = qi*ones(1,length(l));
pl = pi - (lambda*qi^2)/(2*rho*D*A^2).*l;

x1 = 1/2*sqrt(L*C)*pl+1/2*L*ql;
x2 = -1/2*sqrt(L*C)*pl+1/2*L*ql;


% w oparciu o modele transmitancyjne

omega = 0.0000000001;             % stan ustalony


% uk³ad zlinearyzowany - silnie sprzê¿ony 

[dum,dum,dum,dum,G11wi,G12wi,G21wi,G22wi] = freq_resp_trans_line(omega,l,R,L,G,C);

pl_t = G11wi * pi + G12wi * qi;
ql_t = G21wi * pi + G22wi * qi;

% uk³ad rozprzê¿ony 

E = [C 0;0 L];
F = [0 1;1 0];
H = [0 0;0 -R];

% diagonalizacja
[S,LAMBDA]=eig(F*inv(E));
K = inv(S)*H*inv(E)*S;

[dum,dum,dum,dum,G11i,G12i,G21i,G22i] = freq_resp_decoupled(omega,l,LAMBDA,K);

% w celu wyznaczenia transformowanych wejœæ brzegowych

SE = inv(S)*E;
ES = inv(E)*S;
SEd = diag(diag(SE)); %diagonal of SE
SEa = fliplr(diag(diag(fliplr(SE)))); %antidiagonal of SE
ESd = diag(diag(ES)); %diagonal of ES
ESa = fliplr(diag(diag(fliplr(ES)))); %antidiagonal of ES

% opposite boundary matrix
Go = [ G11i(end) G12i(end); G21i(1) G22i(1)];
F = SEa*(ESd*Go + ESa);
I = eye(2);
uw = [pi; qi];         % wektor oryginalnych wejœæ
xi = inv(I-F)*SEd*uw;  % wektor transformowanych wejœæ

G11wi_ = zeros(1,length(l)) ;
G12wi_ = zeros(1,length(l)) ;
G21wi_ = zeros(1,length(l)) ;
G22wi_ = zeros(1,length(l)) ;

x1i = xi(1);
x2i = xi(2);

x1_t = G11i * x1i + G12i * x2i;
x2_t = G21i * x1i + G22i * x2i;

figure(1)
[AX,H1,H2] = plotyy(l,pl,l,ql);
set(H1,'Color','k','LineStyle','-')
set(H2,'Color','k','LineStyle','-')
set(AX(1),'YColor','k')
set(AX(1),'YLim',[4.00E5 5.05E5])
set(AX(1),'YTick',[4.00E5:0.1E5:5.05E5])
set(AX(2),'YColor','k')
set(AX(2),'YLim',[26 29])
set(AX(2),'YTick',[20:0.5:29])

hold on
[AX,H1,H2] = plotyy(l,pl_t,l,ql_t);
set(H1,'Color','k','LineStyle','--')
set(H2,'Color','k','LineStyle','--')
set(AX(1),'YColor','k')
set(AX(1),'YLim',[4.00E5 5.05E5])
set(AX(1),'YTick',[4.00E5:0.1E5:5.05E5])
set(AX(2),'YColor','k')
set(AX(2),'YLim',[26 29])
set(AX(2),'YTick',[26:0.5:29])

set(get(AX(1),'Ylabel'),'String','p(l) [Pa]') 
set(get(AX(2),'Ylabel'),'String','q(l) [kg/s]') 

xlabel('l [m]')
text(100,4.9E5,'p(l)')
text(100,4.6E5,'q(l)')
legend('exact solution','linear model')

% zmienne transformowane
figure(2)
plot(l,x1,'k-')
hold on
plot(l,x1_t,'k--')
plot(l,x2,'k-')
plot(l,x2_t,'k--')
legend('exact solution','linear model')

xlabel('l [m]')
ylabel('x_1(l),x_2(l) [kg/(m^2s)]')
text(100,580,'x_1(l)')
text(100,300,'x_2(l)')


