
function [CM] = COM(Mass,p)
%COM Summary of this function goes here
%   Calculate the centre of mass
%input arguments:
%   Mass    : (1xN) mass vector
%   p       : (3xN) position vector
%output arguments:
%   CM      : (double) 
    CM = nansum(repmat(Mass, [3 1]).*p,2)/sum(Mass);
end

