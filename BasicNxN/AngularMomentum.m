function [L] = AngularMomentum(p,N,Mass,v)
%ANGULARMOMENTUM Summary of this function goes here
%   Detailed explanation goes here
    CM = COM(Mass,p);
    Arm = p - repmat(CM,[1,N]);
    
    %make a momentum vector
    m = v.*repmat(Mass,[3,1]);
    L(1) = nansum(Arm(2,:).*m(3,:) - Arm(3,:).*m(2,:));
    L(2) = nansum(Arm(3,:).*m(1,:) - Arm(1,:).*m(3,:));
    L(3) = nansum(Arm(1,:).*m(2,:) - Arm(2,:).*m(1,:));
    
    %return the [Lx Ly Lz] vector.
end

