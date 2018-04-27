function [D,R] = dispVec(p,N)
%DISPVEC Summary of this function goes here
%   Detailed explanation goes here
    p = permute(p,[3,2,1]); %Make the 'xyz' the third dimension. 
    %The first dimension singleton and the 2nd dimension the particles
    repVec = repmat(p,N,1,1); %Repeat the vectors;
    
    D = permute(repVec,[2,1,3])-repVec; %Transpose and subtract the vectors
    %That way subtracting all the different combinations
    
    %Calculate the range between each particle (stored in NxN matrix).
    R = sqrt(D(:,:,1).^2 + D(:,:,2).^2 + D(:,:,3).^2);
end
