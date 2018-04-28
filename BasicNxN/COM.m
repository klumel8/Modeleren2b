
function [CM] = COM(Mass,p)
%COM Summary of this function goes here
%   Calculate the centre of mass
    CM = nansum(repmat(Mass, [3 1]).*p,2)/sum(Mass);
end

