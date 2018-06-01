function [D,R] = dispVec(p,N)
%DISPVEC calculate the distance(cartesian coordinates) and the
%distance(norm) between the particles
%input arguments:
%   p       : (3xN) position vector
%   N       : (scalar) number of particles
%output arguments:
%   D       : (NxNx3) distance: (i,j,k) = distance between the ith and jth
%               particle, the kth cartesian coordinate(x,y,z)
%   R       : (NxN) distance:(i,j) = distance between the ith and jth
%               particle
%no nonstandard functions required
    p = permute(p,[3,2,1]); %Make the 'xyz' the third dimension.  (1xNx3)
    %The first dimension singleton and the 2nd dimension the particles
    repVec = repmat(p,N,1,1); %Repeat the vectors; (NxNx3)

    D = permute(repVec,[2,1,3])-repVec; %Transpose and subtract the vectors
    %That way subtracting all the different combinations(NxNx3)
    
    %Calculate the range between each particle (stored in NxN matrix).
    R = realsqrt(sum(D.^2,3));
end

