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
text(0.5,96,'x_1(l)','FontSize',14)
text(0.5,103,'x_{1,n}','FontSize',14)
if vs>0
    %text(0.5,50,'\vartheta_s(l)','FontSize',14)
    text(0.5,56,'x_2(l)','FontSize',14)
    text(0.5,48,'x_{2,n}','FontSize',14)
    %text(0.5,70,'\vartheta_w(l)','FontSize',14)
else
    %text(0.5,71,'\vartheta_s(l)','FontSize',14)
    text(0.5,71,'x_2(l)','FontSize',14)
    text(0.5,64,'x_{2,n}','FontSize',14)
    %text(0.5,83,'\vartheta_w(l)','FontSize',14)
end
axis([0 5 43 105])

% approximation models with different numbers of sections

set(0, 'DefaultLineLineWidth', 0.8);
set(0, 'DefaultAxesLineStyleOrder','--|--|--|--')

i = 2;

for N = [2 11 101 1001]
  switch N    
    case 2
        s = [0 0 1];
        p = '*';
    case 11
        s = [1 0 0];
        p = 'o';
    case 101
        s = [0 0.5 0];
        p = 'x';
    case 1001
        s = [1 0 1]
        p = '.';
        N = 101;
  end
  
  ln = [0:L/N:L];      % calculating steady-state profiles
  [GttNw,GtsNw,GstNw,GssNw] = calc_exchanger_approx_freq_resp(omega,ln,v,k); 
  
  theta_t = GttNw*theta_t_i + GtsNw*theta_s_i;  
  theta_s = GstNw*theta_t_i + GssNw*theta_s_i;    

  % vt>0, vs<0

   ln = [0:L/N:L-L/N];            % drawing theta_t profile
   plot([L-L/N L],[theta_t(end) theta_t(end)],'Color',s)
   stairs(ln,[theta_t_i theta_t],'Marker',p,'Color',s,'LineWidth',0.7)
   
   ln = [L:-L/N:L/N];             % drawing theta_s profile
   plot([0 L/N],[theta_s(1) theta_s(1)],'Color',s)
   h = stairs(ln,[theta_s_i theta_s(end:-1:1)],'Marker',p,'Color',s,'LineWidth',0.7)
  
  % vt>0, vs>0

%    ln = [0:L/N:L-L/N];            % drawing theta_t profile
%    plot([L-L/N L],[theta_t(end) theta_t(end)],'Color',s)
%    stairs(ln,[theta_t_i theta_t],'Marker',p,'Color',s,'LineWidth',0.7)
%  
%    ln = [0:L/N:L-L/N];            % drawing theta_s profile
%    plot([L-L/N L],[theta_s(end) theta_s(end)],'Color',s)
%    h = stairs(ln,[theta_s_i theta_s],'Marker',p,'Color',s,'LineWidth',0.7)  

   %ln = [L/N:L/N:L-L/N];          % drawing theta_w profile
   %theta_w = (k2*theta_t + k2*theta_s)/(k2+k3);
   %h = stairs(ln,theta_w,'Marker',p,'Color',s,'LineWidth',0.7)
  H(i) = h;
  i = i+1
end

legend(H,'PDE','N=1','N=10','N=100','N=1000','Location','northoutside','Orientation','horizontal')
