function [e a] = eccentricity_sma(p,v,m)
%Calculate the center of mass of the given particlesx
%input arguments:
%   p       : (3xN) position vector
%   v       : (3xN) velocity vector
%output arguments:
%   e       : (1xN) eccentricity vector
%no nonstandard functions required
    G = 6.67408*10^-11; % [Nm^2kg^-2]

    mu = G * m;
    e = (vecnorm(v).^2 ./ mu - vecnorm(p).^(-1)) * p - p .* v .* v ./ mu;
    
    eps = vecnorm(v).^2 / 2 - mu ./ vecnorm(r);
    
    a = mu ./ (2 * eps);
end