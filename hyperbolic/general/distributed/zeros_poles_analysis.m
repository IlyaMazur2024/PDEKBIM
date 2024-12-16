clear all
close all

L = 5;          
                
lambda1 = 1;    
lambda2 = 0.2;  
                
k11 = -0.05;     
k12 =  0.05;
k21 =  0.05;
k22 = -0.05;

real_line = [-0.1:0.01:0.15];
complex_line = [-0.8:0.01:0.8];

re_s = repmat(real_line,length(complex_line),1);
im_s = repmat(complex_line',1,length(real_line));

p11 = (k11-(re_s+i*im_s))/lambda1;
p12 = k12/lambda1*ones(length(complex_line),length(real_line));
p21 = k21/lambda2*ones(length(complex_line),length(real_line));
p22 = (k22-(re_s+i*im_s))/lambda2;

alpha = 1/2*(p11+p22);
beta = 1/2*sqrt((p11-p22).^2+4*p12.*p21);

% phi1, phi22 - wartoœci w³asne macierzy P

phi1 = alpha + beta;
phi2 = alpha - beta;

figure(1)
mesh(re_s,im_s,abs(phi1-phi2))
xlabel('Re(s)')
ylabel('Im(s)')
zlabel('|\phi_1(s)-\phi_2(s)|')
axis([-0.1 0.1 -1 1 0 2])
min(min(abs(phi1-phi2)))

figure(2)
mesh(re_s,im_s,abs(exp(phi1*L)-exp(phi2*L)))
xlabel('Re(s)')
ylabel('Im(s)')
zlabel('|exp(\phi_{ 1}(s)l)-exp(\phi_{ 2}(s)l)|')
min(min(abs(exp(phi1*L)-exp(phi2*L))))

figure(3)
mesh(re_s,im_s,abs(exp(phi2*L).*(phi1-p22)-exp(phi1*L).*(phi2-p22)))
xlabel('Re(s)')
ylabel('Im(s)')
min(min(abs(exp(phi2*L).*(phi1-p22)-exp(phi1*L).*(phi2-p22))))

figure(4)
mesh(re_s,im_s,abs(exp(phi2*L).*(phi1-p22)-exp(phi1*L).*(phi2-p22)))
xlabel('Re(s)')
ylabel('Im(s)')
min(min(abs(exp(phi2*L).*(phi1-p22)-exp(phi1*L).*(phi2-p22))))


% pole/zero phi1=phi2
p1=(k22*lambda1-k11*lambda2+2*sqrt(-k12*k21*lambda1*lambda2))/(lambda1-lambda2)
p2=(k22*lambda1-k11*lambda2+2*sqrt(-k12*k21*lambda1*lambda2))/(lambda1-lambda2)


% zeros exp(phi1)=exp(phi2)
k=2

z1=(k22*lambda1-k11*lambda2)/(lambda1-lambda2)+...
    2*sqrt(-lambda1^2*lambda2^2*k^2*pi^2-lambda1*lambda2*k12*k21*L^2)/((lambda1-lambda2)*L)

z2=(k22*lambda1-k11*lambda2)/(lambda1-lambda2)-...
    2*sqrt(-lambda1^2*lambda2^2*k^2*pi^2-lambda1*lambda2*k12*k21*L^2)/((lambda1-lambda2)*L)


