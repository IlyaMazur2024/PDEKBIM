%function [ SAB11,SAB12,SAB21,SAB22 ] = redhefferstar( SA11,SA21,SA12,SA22,SB11,SB21,SB12,SB22)
%function [g11,g12,g21,g22] = redhefferstar( g11A,g12A,g21A,g22A,g11B,g12B,g21B,g22B)
function [G] = redhefferstar( g11A,g12A,g21A,g22A,g11B,g12B,g21B,g22B)

% redhefferstar product combines two scattering matrices to form a overall
% scattering matrix. It is used in forming scattering matrices of
% dielectric stacks in transfer matrix method
% SA and SB are  scattering matrices of two different layers
% and this function outputs
% SAB which is the combined scaterring matrix of two layers

N=length(g11A);
I=eye(N);

g11 = g11B *(I-g12A*g21B)^(-1)*g11A;
g12 = g11B*(I-g12A*g21B)^(-1)*g12A*g22B+g12B;
g21 = g21A + g22A*(I-g21B*g12A)^(-1)*g21B*g11A;
g22 = g22A*(I-g21B*g12A)^(-1)*g22B;

G = [g11 g12; g21 g22];

end

