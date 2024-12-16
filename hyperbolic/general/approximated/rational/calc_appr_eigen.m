
close all
clear all
init_hyperb
set(0,'DefaultLineLineWidth',1.5)
set(0,'defaultAxesFontSize',14)

N=round(logspace(0,3,200));
%N=[1:10:1000];
ev = zeros(2,length(N)); 
i=1;

for n = N
    dl = L/n;
    An = [-lambda1/dl+k11 k12; k21 -lambda2/dl+k22];
    ev(:,i) = eig(An);
    i=i+1;
end


figure(1)
loglog(N,ev(1,:),'k.')
%plot(N,ev(1,:),'k.')
hold on
loglog(N,ev(2,:),'r.')
%plot(N,ev(2,:),'r.')
grid on
xlabel('N')
%ylabel('log_{10}(eig\{A_n\})')
ylabel('\sigma_{(1,2),n}')
legend('\sigma_{1,n}','\sigma_{2,n}')