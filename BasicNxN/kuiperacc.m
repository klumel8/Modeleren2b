function a = kuiperacc(p, p_k, Mass)
% Calculate the acceleration of all the particles based on the
% another system of particles with mass
%input arguments:
%   p       : (3xN) position matrix
%   Mass    : (1xN) mass vector
%   G       : (scalar) gravitational constant
%   N       : (scalar)number of particles
%output arguments:
%   a       : (3xN) acceleration matrix
%required functions(non-standard):
%   dispVec

    G = 6.67408*10^-11; % [Nm^2kg^-2]
    N = size(p,2);
    N_k = size(p_k,2);
    a = zeros(3,N_k,N);
    %get the distances between the particles
    for i = 1
    D = repmat(p(:,i), [1,N_k]) - p_k;
    R = vecnorm(D);

    %Make a simple mathematical expression for everything distance related,
    %later to be used in the force formulae.
    temp= repmat(R, [3,1]);
    %gpuTemp = gpuArray(temp);
    %gpuD = gpuArray(D);
    %gpuDisp = gpuD./(gpuTemp.*gpuTemp.*gpuTemp);
    
    dispComp = D./(temp.^3);
    
    %calculate the force in each direction (x,y,z).
    %no minus sign because somewhere I misplaced a minus sign (dont worry
    %I anticipated this)
    a(:,:,i) = G*repmat(Mass(i), [3,1]).*dispComp;
    end
    a = sum(a,3);
end

