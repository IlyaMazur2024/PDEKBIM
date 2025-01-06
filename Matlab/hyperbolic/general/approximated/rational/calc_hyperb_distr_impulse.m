
function [g11,g12,g21,g22] = calc_hyper_distr_impulse(LAMBDA,K,L,t,l,N)
    
% =========================================================================
% 
% Function m-file name: analytical_impulse_responses.m
%
% Purpose: function calculates spatiotemporal impulse reponses of 
% the following weakly coupled hyperbolic system: 
%
% dx1(l,t)/dt + lambda1*dx1(l,t)/dl = k11*x1(l,t) + k12*x2(l,t)		       
% dx2(l,t)/dt + lambda2*dx2(l,t)/dl = k21*x1(l,t) + k22*x2(l,t)		       
%
% assuming the collocated boundary inputs: 
% x1(0,t), x2(0,t) for lambda1>0, lambda2>0, 
% or the anti-collocated ones: 
% x1(0,t), x2(L,t) for lambda1>0, lambda2>0, 
% as well as zero initial conditions: x1(l,0)=0, x2(l,0)=0.
%
% The impulse responses are obtained from the inverse Laplace transform 
% of the transfer functions of the system.
% 
% Syntax: [g11,g12,g21,g22]=calc_hyper_distr_impulse(LAMBDA,K,L,t,l,N)
%
% Function parameters:
% LAMBDA - 2 x 2 diagonal matrix of eigenvalues, 
% LAMBDA = diag(lambda1,lambda2) where lambda1>lambda2
% K - 2 x 2 coupling matrix, 
% t - vector of time instants,
% l - vector of spatial positions,
% L - length of the spatial domain.
% N - number of waves included in the impulse response 
%     (matters only for incongruent inputs)
%
% Output values:
% g11, g12, g21, g22 - Nt x Nl matrices containing spatiotemporal impulse
%                      responses of the corresponding input-output channels, 
%                      where: Nt = lenght(t), Nl = length(l).
%
% =========================================================================

lambda1=LAMBDA(1,1);
lambda2=LAMBDA(2,2);

k11 = K(1,1);
k12 = K(1,2);
k21 = K(2,1);
k22 = K(2,2);

l_ = repmat(l,length(t),1);
t_ = repmat(t',1,length(l));

eta = 4*k12*k21*lambda1*lambda2;  

if lambda1>0 && lambda2>0 % congruent boundary conditions  

epsilon11 = -sqrt(eta)/(2*(lambda1-lambda2));
epsilon12 = k12*lambda2/(lambda1-lambda2);
epsilon21 = k21*lambda1/(lambda1-lambda2);
epsilon22 = epsilon11;

tau1_ = l_/lambda1;
tau2_ = l_/lambda2;

kappa1_ = exp(k11*l_/lambda1);
kappa2_ = exp(k22*l_/lambda2);

psi11_1_ = exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau1_)).*...
          sqrt((t_-tau2_)./(t_-tau1_)).*besselj(1,sqrt(eta)/(lambda1-lambda2)*sqrt((t_-tau1_).*(t_-tau2_)));
psi11_2_ = exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau2_)).*...
          sqrt((t_-tau2_)./(t_-tau1_)).*besselj(1,sqrt(eta)/(lambda1-lambda2)*sqrt((t_-tau1_).*(t_-tau2_)));
      
psi12_1_ = exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau1_)).*...
          besselj(0,sqrt(eta)/(lambda1-lambda2).*sqrt((t_-tau1_).*(t_-tau2_)));
psi12_2_ = exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau2_)).*...
          besselj(0,sqrt(eta)/(lambda1-lambda2).*sqrt((t_-tau1_).*(t_-tau2_)));      
      
psi21_1_ = exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau1_)).*...
          besselj(0,sqrt(eta)/(lambda1-lambda2).*sqrt((t_-tau1_).*(t_-tau2_)));
psi21_2_ = exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau2_)).*...
          besselj(0,sqrt(eta)/(lambda1-lambda2).*sqrt((t_-tau1_).*(t_-tau2_)));            
      
psi22_1_ = exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau1_)).*...
          sqrt((t_-tau1_)./(t_-tau2_)).*besselj(1,sqrt(eta)/(lambda1-lambda2).*sqrt((t_-tau1_).*(t_-tau2_)));
psi22_2_ = exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau2_)).*...
          sqrt((t_-tau1_)./(t_-tau2_)).*besselj(1,sqrt(eta)/(lambda1-lambda2).*sqrt((t_-tau1_).*(t_-tau2_)));
      
g11 = epsilon11*(heaviside(t_-tau1_).*kappa1_.*psi11_1_ - heaviside(t_-tau2_).*kappa2_.*psi11_2_);
g12 = epsilon12*(heaviside(t_-tau1_).*kappa1_.*psi12_1_ - heaviside(t_-tau2_).*kappa2_.*psi12_2_);
g21 = epsilon21*(heaviside(t_-tau1_).*kappa1_.*psi21_1_ - heaviside(t_-tau2_).*kappa2_.*psi21_2_);
g22 = epsilon22*(heaviside(t_-tau1_).*kappa1_.*psi22_1_ - heaviside(t_-tau2_).*kappa2_.*psi22_2_);

elseif lambda1>0 && lambda2<0 % incongruent boundary conditions
   
epsilon11 = 1;
epsilon12 = -2/sqrt(-eta)*k12*lambda2;
epsilon21 =  2/sqrt(-eta)*k21*lambda1;
epsilon22 = 1;

g11 = zeros(length(t),length(l));
g12 = zeros(length(t),length(l));
g21 = zeros(length(t),length(l));
g22 = zeros(length(t),length(l));

  for k=0:N  % N - number of waves included in the response
    
    tau1_1_ = (k*L +l_)/lambda1 - k*L/lambda2;
    tau1_2_ = (k+1)*L/lambda1 - ((k+1)*L-l_)/lambda2;

    kappa1_1_ = exp(k11*(k*L +l_)/lambda1 - k22*k*L/lambda2);
    kappa1_2_ = exp(k11*(k+1)*L/lambda1 - k22*((k+1)*L-l_)/lambda2);
    
    mu1_1_ = -(2*k*L+l_)/(2*lambda1*lambda2);
    mu1_2_ = -(2*(k+1)*L-l_)/(2*lambda1*lambda2);
    
    tau2_1_ = k*L/lambda1 - ((k+1)*L-l_)/lambda2;
    tau2_2_ = (k*L+l_)/lambda1 - (k+1)*L/lambda2;
    
    kappa2_1_ = exp(k11*k*L/lambda1 - k22*((k+1)*L-l_)/lambda2);
    kappa2_2_ = exp(k11*(k*L+l_)/lambda1 - k22*(k+1)*L/lambda2);
    
    mu2_1_ = -((2*k+1)*L-l_)/(2*lambda1*lambda2);
    mu2_2_ = -((2*k+1)*L+l_)/(2*lambda1*lambda2);
    
    psi_11_1_ = 1./(t_-tau1_1_).*exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau1_1_)).*...
                (sqrt(-eta)*mu1_1_.*((t_-tau1_1_)./(t_-tau1_1_+2*mu1_1_*(lambda1-lambda2))).^(k+1/2).*...
                besseli(2*k+1,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau1_1_).*((t_-tau1_1_)+2*mu1_1_*(lambda1-lambda2))))+...
                2*k*((t_-tau1_1_)./(t_-tau1_1_+2*mu1_1_*(lambda1-lambda2))).^k.*...
                besseli(2*k,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau1_1_).*(t_-tau1_1_+2*mu1_1_*(lambda1-lambda2)))));
           
    psi_11_2_ = 1./(t_-tau1_2_).*exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau1_2_)).*...
                (sqrt(-eta)*mu1_2_.*((t_-tau1_2_)./(t_-tau1_2_+2*mu1_2_*(lambda1-lambda2))).^(k+3/2).*...
                besseli(2*k+3,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau1_2_).*((t_-tau1_2_)+2*mu1_2_*(lambda1-lambda2))))+...
                (2*k+2)*((t_-tau1_2_)./(t_-tau1_2_+2*mu1_2_*(lambda1-lambda2))).^(k+1).*...
                besseli(2*k+2,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau1_2_).*(t_-tau1_2_+2*mu1_2_*(lambda1-lambda2)))));

    psi_12_1_ = 1./(t_-tau2_1_).*exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau2_1_)).*...
                (sqrt(-eta)*mu2_1_.*((t_-tau2_1_)./(t_-tau2_1_+2*mu2_1_*(lambda1-lambda2))).^(k+1).*...
                besseli(2*k+2,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau2_1_).*((t_-tau2_1_)+2*mu2_1_*(lambda1-lambda2))))+...
                (2*k+1)*((t_-tau2_1_)./(t_-tau2_1_+2*mu2_1_*(lambda1-lambda2))).^(k+1/2).*...
                besseli(2*k+1,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau2_1_).*(t_-tau2_1_+2*mu2_1_*(lambda1-lambda2)))));            

    psi_12_2_ = 1./(t_-tau2_2_).*exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau2_2_)).*...
                (sqrt(-eta)*mu2_2_.*((t_-tau2_2_)./(t_-tau2_2_+2*mu2_2_*(lambda1-lambda2))).^(k+1).*...
                besseli(2*k+2,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau2_2_).*((t_-tau2_2_)+2*mu2_2_*(lambda1-lambda2))))+...
                (2*k+1)*((t_-tau2_2_)./(t_-tau2_2_+2*mu2_2_*(lambda1-lambda2))).^(k+1/2).*...
                besseli(2*k+1,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau2_2_).*(t_-tau2_2_+2*mu2_2_*(lambda1-lambda2)))));   

    psi_21_1_ = 1./(t_-tau1_1_).*exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau1_1_)).*...
                (sqrt(-eta)*mu1_1_.*((t_-tau1_1_)./(t_-tau1_1_+2*mu1_1_*(lambda1-lambda2))).^(k+1).*...
                besseli(2*k+2,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau1_1_).*((t_-tau1_1_)+2*mu1_1_*(lambda1-lambda2))))+...
                (2*k+1)*((t_-tau1_1_)./(t_-tau1_1_+2*mu1_1_*(lambda1-lambda2))).^(k+1/2).*...
                besseli(2*k+1,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau1_1_).*(t_-tau1_1_+2*mu1_1_*(lambda1-lambda2)))));                       

    psi_21_2_ = 1./(t_-tau1_2_).*exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau1_2_)).*...
                (sqrt(-eta)*mu1_2_.*((t_-tau1_2_)./(t_-tau1_2_+2*mu1_2_*(lambda1-lambda2))).^(k+1).*...
                besseli(2*k+2,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau1_2_).*((t_-tau1_2_)+2*mu1_2_*(lambda1-lambda2))))+...
                (2*k+1)*((t_-tau1_2_)./(t_-tau1_2_+2*mu1_2_*(lambda1-lambda2))).^(k+1/2).*...
                besseli(2*k+1,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau1_2_).*(t_-tau1_2_+2*mu1_2_*(lambda1-lambda2)))));                                    
            
    psi_22_1_ = 1./(t_-tau2_1_).*exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau2_1_)).*...
                (sqrt(-eta)*mu2_1_.*((t_-tau2_1_)./(t_-tau2_1_+2*mu2_1_*(lambda1-lambda2))).^(k+1/2).*...
                besseli(2*k+1,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau2_1_).*((t_-tau2_1_)+2*mu2_1_*(lambda1-lambda2))))+...
                2*k*((t_-tau2_1_)./(t_-tau2_1_+2*mu2_1_*(lambda1-lambda2))).^k.*...
                besseli(2*k,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau2_1_).*(t_-tau2_1_+2*mu2_1_*(lambda1-lambda2)))));
            
    psi_22_2_ = 1./(t_-tau2_2_).*exp(-(k11*lambda2-k22*lambda1)/(lambda1-lambda2)*(t_-tau2_2_)).*...
                (sqrt(-eta)*mu2_2_.*((t_-tau2_2_)./(t_-tau2_2_+2*mu2_2_*(lambda1-lambda2))).^(k+3/2).*...
                besseli(2*k+3,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau2_2_).*((t_-tau2_2_)+2*mu2_2_*(lambda1-lambda2))))+...
                (2*k+2)*((t_-tau2_2_)./(t_-tau2_2_+2*mu2_2_*(lambda1-lambda2))).^(k+1).*...
                besseli(2*k+2,sqrt(-eta)/(lambda1-lambda2)*sqrt((t_-tau2_2_).*(t_-tau2_2_+2*mu2_2_*(lambda1-lambda2)))));           
            
    g11 = g11 + epsilon11*(heaviside(t_-tau1_1_).*kappa1_1_.*psi_11_1_ - heaviside(t_-tau1_2_).*kappa1_2_.*psi_11_2_);
    g12 = g12 + epsilon12*(heaviside(t_-tau2_1_).*kappa2_1_.*psi_12_1_ - heaviside(t_-tau2_2_).*kappa2_2_.*psi_12_2_);
    g21 = g21 + epsilon21*(heaviside(t_-tau1_1_).*kappa1_1_.*psi_21_1_ - heaviside(t_-tau1_2_).*kappa1_2_.*psi_21_2_);
    g22 = g22 + epsilon22*(heaviside(t_-tau2_1_).*kappa2_1_.*psi_22_1_ - heaviside(t_-tau2_2_).*kappa2_2_.*psi_22_2_);
    
  end 
end







