function [e, a] = eccentricity_sma(p,v,m)
%Calculate the eccentricity and the semi major axis
%input arguments:
%   p       : (3xN) position vector
%   v       : (3xN) velocity vector
%   m       : (1xN) Mass vector
%output arguments:
%   e       : (1xN) eccentricity
%   a       : (1xN) semi major axis
%no nonstandard functions required
    G = 6.67408*10^-11; % [Nm^2kg^-2]

    mu = G * m(1);
    eps = vecnorm(v).^2 / 2 - mu./vecnorm(p); 
    a = -mu ./ (2 * eps);
    e = sqrt(1 - sum(cross(p,v).^2)./(mu*a));
end