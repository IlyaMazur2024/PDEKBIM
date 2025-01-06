clear
syms L c R G positive 
syms E F W real
syms lambda1 lambda2 k11 k12 k21 k22 real
syms s1 s2 s b c alpha beta q1 q2 r11 r22 complex
syms r12 r21 l LL real


% obliczanie transmitanacji operatorowej linii elektrycznej
% w oparciu o model rozsprzê¿ony

E = [c 0;0 L];
%E = [L 0;0 c];
%E = [0 c;L 0];

F = [0 1;1 0];
%F = [1 0;0 1];

%W = [0 0;0 0];
W = [-R 0;0 -G];
%W = [0 0;0 -R];   


% diagonalizacja

[S,LAMBDA]=eig(F*inv(E));
LAMBDA=rot90(rot90(LAMBDA))
S=fliplr(S)
%lambda1 =  1/sqrt(L*c);
%lambda2 = -1/sqrt(L*c);

K = inv(S)*W*inv(E)*S;
%k11 = -G/(2*c)-R/(2*L); 
%k12 = -G/(2*c)+R/(2*L); 
%k21 = -G/(2*c)+R/(2*L); 
%k22 = -G/(2*c)-R/(2*L); 

I = eye(2);
P = inv(LAMBDA)*(K-s*I);
p11 = P(1,1);
p12 = P(1,2);
p21 = P(2,1);
p22 = P(2,2);

[PSI,PHI] = eig(P);
p=diag(PHI);
phi1=p(1);
phi2=p(2);


% =========================================================================

% transmitancje operatorowe w przypadku obydwu warunków
% brzegowych okreœlonych dla pocz¹tku linii

%l = 1;

% transmitancje operatorowe (widmowe) uk³adu rozsprzê¿onego 
% z wykorzystaniem funkcji wyk³adniczych

A11 = (phi1-p22)./(phi1-phi2);
B11 = (phi2-p22)./(phi1-phi2);
A12 = p12./(phi1-phi2);
A21 = p21./(phi1-phi2);
A22 = (phi1-p11)./(phi1-phi2);
B22 = (phi2-p11)./(phi1-phi2);

Gx11 = A11 .* exp(phi1*l) - B11 .* exp(phi2*l);  
Gx12 = A12.*(exp(phi1.*l) - exp(phi2.*l));       
Gx21 = A21.*(exp(phi1.*l) - exp(phi2.*l));       
Gx22 = A22.*exp(phi1*l) - B22.* exp(phi2*l);    

Gx = [Gx11 Gx12;Gx21 Gx22];
% konwersja do transmitancji operatorowych 
% uk³adu oryginalnego (sprzê¿onego):
Gw = inv(E)*S*Gx*inv(S)*E;
pretty(simple(Gw))

% =========================================================================

% transmitancje operatorowe w przypadku obydwu warunków
% brzegowych przeciwstawnych
A11 = (exp(phi2*LL).*(phi1-p22)) ./ (exp(phi2*LL).*(phi1-p22) - exp(phi1*LL).*(phi2-p22));
B11 = (exp(phi1*LL).*(phi2-p22)) ./ (exp(phi2*LL).*(phi1-p22) - exp(phi1*LL).*(phi2-p22));
A12 = p12 ./ (exp(phi2*LL).*(phi2-p11) - exp(phi1*LL).*(phi1-p11));
A21 = p21*exp(phi2*LL) ./ (exp(phi2*LL).*(phi1-p22) - exp(phi1*LL).*(phi2-p22));
B21 = p21*exp(phi1*LL) ./ (exp(phi2*LL).*(phi1-p22) - exp(phi1*LL).*(phi2-p22));
A22 = (phi2-p11) ./ (exp(phi2*LL).*(phi2-p11) - exp(phi1*LL).*(phi1-p11));
B22 = (phi1-p11) ./ (exp(phi2*LL).*(phi2-p11) - exp(phi1*LL).*(phi1-p11));

Gx11_ = A11 .* exp(phi1*l) - B11 .* exp(phi2*l);  
Gx12_ = A12 .* (exp(phi2.*l) - exp(phi1.*l));       
Gx21_ = A21 .* exp(phi1.*l) - B21 .* exp(phi2.*l) ;       
Gx22_ = A22 .* exp(phi2*l) - B22.* exp(phi1*l);    

Gx_ = [Gx11_ Gx12_;Gx21_ Gx22_];
% konwersja do transmitancji operatorowych 
% uk³adu oryginalnego (sprzê¿onego):

SE = inv(S)*E;
ES = inv(E)*S;

SEd = diag(diag(SE)); %diagonal of SE
SEa = fliplr(diag(diag(fliplr(SE)))); %antidiagonal of SE
ESd = diag(diag(ES)); %diagonal of ES
ESa = fliplr(diag(diag(fliplr(ES)))); %antidiagonal of ES

G11b_ = A11 .* exp(phi1*LL) - B11 .* exp(phi2*LL);  
G12b_ = A12.*(exp(phi2.*LL) - exp(phi1.*LL));       
G21b_ = A21.* exp(phi1.*0) - B21 .* exp(phi2.*0) ;       
G22b_ = A22.*exp(phi2*0) - B22.* exp(phi1*0);    

Gb_ = [G11b_ G12b_; G21b_ G22b_];
I = eye(2);

%Gw_ = ES*Gx_*inv(I-SEa*(ESa+ESd*Gb_))*SEd;
%pretty(simple(Gw_))

% nowy krótszy wzór
Gw_ = ES*Gx_*inv(ESd+ESa*Gb_);
%pretty(simple(Gw_))

