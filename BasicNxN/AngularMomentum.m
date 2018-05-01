function [L] = AngularMomentum(p,N,Mass,v)
%ANGULARMOMENTUM calculate the total angular momentum of the particles
%input arguments:
%   p       : (3xN) position vector
%   N       : (scalar) number of particles
%   Mass    : (1xN) mass vector
%   v       : (3xN) velocity vector
%output arguments:
%   L       : (1x3) angular momentum vector
%required functions (non-standard):
%   COM
    %get the center of mass and calculate the relative distance to the each
    %particle from the center of mass
    CM = COM(Mass,p);
    Arm = p - repmat(CM,[1,N]);
    
    %make a momentum vector
    m = v.*repmat(Mass,[3,1]);
    
    %calculate the angular momentum
    L(1) = nansum(Arm(2,:).*m(3,:) - Arm(3,:).*m(2,:));
    L(2) = nansum(Arm(3,:).*m(1,:) - Arm(1,:).*m(3,:));
    L(3) = nansum(Arm(1,:).*m(2,:) - Arm(2,:).*m(1,:));
    
    %return the [Lx Ly Lz] vector.
end

