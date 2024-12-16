% ===========================================================================
%
% Comaprison of the steady-state temperature distribution calculated from 
% the PDE model of the parallel-flow heat exchanger and the approximate models 
% with different numbers of spatial sections
close all
clear all

init_exchanger                 % initialization of the PDE model parameters
omega = 0;                     % steady state conditions
l = [0:0.1:L];                 % entire exchanger length is considered

theta_t_i = 100;         % boundary conditions (constant inlet temperatures)
theta_s_i = 50;           
 
% calculation of the steady-state transfer functions
% from the PDE-based irrational transfer function model
[gtt0,gts0,gst0,gss0] = calc_exchanger_distr_freq_resp(omega,l,v,k);

% exact steady-state profiles 
theta_t_l = gtt0*theta_t_i + gts0*theta_s_i; 
theta_s_l = gst0*theta_t_i + gss0*theta_s_i; 
theta_w_l = (k2*theta_t_l + k2*theta_s_l)/(k2+k3);

H = gobjects(1,4);
set(0,'DefaultLineLineWidth',1)
set(0,'defaultAxesFontSize',14)

figure(1)
plot(l,theta_t_l,'k-')
hold on
%plot(l,theta_w_l,'k-')
h=plot(l,theta_s_l,'k-');
H(1)=h;
xlabel('l','FontSize',14)
%ylabel('\vartheta(l)','FontSize',14)
ylabel('x_1, x_2','FontSize',14)
grid on
%text(0.5,102,'\vartheta_t(l)','FontSize',14)
text(0.5,102,'x_1(l)','FontSize',14)
text(0.5,97,'x_{1,n}(l)','FontSize',14)
if vs>0
    %text(0.5,50,'\vartheta_s(l)','FontSize',14)
    text(0.5,50,'x_2(l)','FontSize',14)
    text(0.5,43,'x_{2,n}(l)','FontSize',14)
    %text(0.5,70,'\vartheta_w(l)','FontSize',14)
else
    %text(0.5,71,'\vartheta_s(l)','FontSize',14)
    text(0.5,71,'x_2(l)','FontSize',14)
    text(0.5,64,'x_{2,n}(l)','FontSize',14)
    %text(0.5,83,'\vartheta_w(l)','FontSize',14)
end
axis([0 5 48 105])

% approximation models with different numbers of sections

set(0, 'DefaultLineLineWidth', 0.8);
set(0, 'DefaultAxesLineStyleOrder','--|--|--|--')

i = 2;
%for N = [2 10 100 1000]
for N = [9 99 999]
  switch N    
    %case 2
    case 2
        s = [0 0 1];
        p = '*';
        load G1w_exchanger_par
        N=1;
    %case 10
    case 9
        s = [1 0 0];
        p = 'o';
        load G10w_exchanger_par
        N=9;
    %case 100
    case 99
        s = [0 0.5 0];
        p = 'x';
        load G100w_exchanger_par
        N=99;
    case 999
        s = [1 0 1]
        p = '.';
        N=99;
  end
  %ln = [0:L/N:L];      % calculating steady-state profiles
  %[GttNw,GtsNw,GstNw,GssNw] = calc_exchanger_approx_freq_resp(omega,ln,v,k); 
  
  GttNw = real(gttw_hat(1,:));
  GtsNw = real(gtsw_hat(1,:));
  GstNw = real(gstw_hat(1,:));
  GssNw = real(gssw_hat(1,:));
  
  theta_t = GttNw*theta_t_i + GtsNw*theta_s_i;  
  theta_s = GstNw*theta_t_i + GssNw*theta_s_i;    

  ln = [0:L/N:L-L/N];            % drawing theta_t profile
  plot([L-L/N L],[theta_t(end) theta_t(end)],'Color',s)
  stairs(ln,[theta_t_i theta_t],'Marker',p,'Color',s,'LineWidth',0.7)
  
  %ln = [L:-L/N:L/N];             % drawing theta_s profile
  plot([0 L/N],[theta_s(1) theta_s(1)],'Color',s)
  h = stairs(ln,[theta_s_i theta_s(end:-1:1)],'Marker',p,'Color',s,'LineWidth',0.7)
     
  
  %ln = [L/N:L/N:L-L/N];          % drawing theta_w profile
  %theta_w = (k2*theta_t + k2*theta_s)/(k2+k3);
  %h = stairs(ln,theta_w,'Marker',p,'Color',s,'LineWidth',0.7)
  H(i) = h;
  i = i+1
end

legend(H,'PDE','N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')
