function [CM] = COM(Mass,p)
%Calculate the center of mass of the given particlesx
%input arguments:
%   Mass    : (1xN) mass vector
%   p       : (3xN) position vector
%output arguments:
%   CM      : (scalar) center of mass
%no nonstandard functions required
    CM = nansum(repmat(Mass, [3 1]).*p,2)/sum(Mass);
end

