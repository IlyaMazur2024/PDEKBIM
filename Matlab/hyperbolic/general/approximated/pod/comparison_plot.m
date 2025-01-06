clc
close all
clear all

O = [1 2 3 5 10 15 20];

EO = [27.93 52.41 65.17 83.32 95.07 98.00 99.05];
CR = [33.83 18.47 12.70 7.82 3.99 2.68 2.01];

O_FULL = [1:0.1:20];
POL_EO = fit(O,EO);
%EO_FULL = ppval(POL_EO,O_FULL);
POL_CR = fit(O,CR);
CR_FULL = ppval(POL_CR,O_FULL);

% O_FULL = [1 2 3 4 5 6 7 8 10 12 14 15 17 20];
% EO_FULL = [27.93 52.41 65.17 77.38 83.32 88.56 91.17 93.20 95.07 96.73 97.64 98.00 98.57 99.05];
% CR_FULL = [33.83 18.47 12.70 9.68 7.82 6.56 5.65 4.96 3.99 3.33 2.86 2.68 2.36 2.01];

C_PCA_EF =   [10.28 8.36 7.15 4.95 2.69 1.71 1.18];
FF_PCA_EF =  [10.28 8.36 7.15 4.95 2.69 1.71 1.18];
GHA_PCA_EF = [10.29 8.37 7.20 5.00 2.78 1.91 1.63];

figure(1)
[AX,H1,H2] = plotyy(O,CR,O,EO);
set(H1,'Color','k','Marker','*','LineStyle','none')
set(H2,'Color','k','Marker','O','LineStyle','none')
legend('CR','E_K')
hold on
[AX,H1,H2] = plotyy(O_FULL,CR_FULL,O_FULL,EO_FULL);
set(H1,'Color','k','Marker','none','LineStyle',':','LineWidth',1)
set(H2,'Color','k','Marker','none','LineStyle',':','LineWidth',1)

set(AX(1),'YColor','k')
set(AX(1),'YLim',[0 35])
set(AX(1),'YTick',[0:5:35])
set(AX(2),'YColor','k')
set(AX(2),'YLim',[20 100])
set(AX(2),'YTick',[20:10:100])

set(get(AX(1),'Ylabel'),'String','CR') 
set(get(AX(2),'Ylabel'),'String','E_K (%)') 
xlabel('O')

%grid on


