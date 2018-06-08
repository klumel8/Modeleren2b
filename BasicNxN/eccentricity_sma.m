function [e, a] = eccentricity_sma(p,v,m)
%Calculate the center of mass of the given particlesx
%input arguments:
%   p       : (3xN) position vector
%   v       : (3xN) velocity vector
%   m       : (1xN) Mass vector
%output arguments:
%   e       : (3xN) eccentricity vector
%   a       : (1xN) semi major axis
%no nonstandard functions required
    G = 6.67408*10^-11; % [Nm^2kg^-2]

    mu = G * m(1);
    e = (sum(v.^2,1) ./ mu - sqrt(sum(p.^2,1)).^(-1)) .* p - dot(p, v,1) .* v ./ mu;
    
    eps = sum(v.^2,1) / 2 - mu ./ sqrt(sum(p.^2,1));
    e = sqrt(e);
    a = -mu ./ (2 * eps);
end