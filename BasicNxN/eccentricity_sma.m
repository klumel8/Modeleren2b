function [e, a] = eccentricity_sma(p,v,m, pp)
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
    p = p - pp(:,1);
    eps = sum(v.^2,1) / 2 - mu./(sum(p.^2,1)).^(0.5); 
    a = -mu ./ (2 * eps);
    e = sqrt(1 - sum(cross(p,v).^2)./(mu*a));
end