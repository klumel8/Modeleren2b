function [a_k,a] = kuiperacc(p, p_k, Mass,Mass_k)
% Calculate the acceleration of all the particles based on the
% another system of particles with mass
%input arguments:
%   p       : (3xN) position matrix
%   Mass    : (1xN) mass vector
%   G       : (scalar) gravitational constant
%   N       : (scalar)number of particles
%output arguments:
%   a_k       : (3xN) acceleration matrix of kuiper particles
%   a         : (3xN) acceleration matrix of planets
%required functions(non-standard):
%   dispVec

    G = 6.67408*10^-11; % [Nm^2kg^-2]


    D_new = permute(p,[1,3,2]) - p_k;%3xN_kxN;
    R_new = sqrt(sum(D_new.^2,1)); %1xN_kxN;


    %Make a simple mathematical expression for everything distance related,
    %later to be used in the force formulae.
%     temp= repmat(R, [3,1]);
    %gpuTemp = gpuArray(temp);
    %gpuD = gpuArray(D);
    %gpuDisp = gpuD./(gpuTemp.*gpuTemp.*gpuTemp);
    
    dispComp_new = D_new./(R_new.^3);
    
    %calculate the force in each direction (x,y,z).
    %no minus sign because somewhere I misplaced a minus sign (dont worry
    %I anticipated this)

    a_k = sum(G*permute(Mass,[1,3,2]).*dispComp_new,3);
    a = sum(-G*Mass_k.*dispComp_new,2);
    a = permute(a,[1,3,2]);


end

