function [p, v] = kuiperbelt(N, p_planet)
%kuiperbelt Creates the initial conditions for the simulation
%
%   Syntax:
%
%       [p v]=kuiperbelt(N)
%
%   Input:
%           * N       = "Number of particles to be simulated" (1)[-]
%                        2 = solar system"                    (1)[-]
%   Output:
%           * p       = "position vector" (3:N)[m]
%           * v       = "velocity vector" (3:N)[m/s]
%
%   Example(s):
%
%       A system of 1.000.000 particles can be creates by
%               [p, v] = kuiperbelt(1e6)
%
%   Dependencies: [-]
%
%   Author: Niels Buijssen
%   Email:  nbuijssen@student.tudelft.nl
%
%   Project: Modelleren 2B Planets
%   Keywords: Planets, Solar system, Initial Conditions
%             
%
%   Version: v1.0 (16-05-2018)
%

%% Change-log
%
%   v1.0 (16-05-2018)
%       * First operational version.
%

%% Input Handling
%   -

%% Main Function

% Define starting parameters 
Mass_sun = 1.988e30; % [kg] mass of sun 
G = 6.67408*10^-11; % [Nm^2kg^-2]
AU = 1.49597871e11;% [m]
% r_low = 42*AU;
% r_high = 48*AU;

% create position and speed vectors
theta = 2*pi*rand(1,N); % create random angles
% r = r_low + (r_high-r_low).*rand(1,N); % create uniformly distributed radii
%r = (2)^(2/3)*4495e9;
 r = ((1)^(2/3))*4495e9;

ecc = rand(1,N)*0.1;

%use defualt gonio functions to make a physically valid semi-major/minor
%axis a,b
a = r;
b = a * sqrt(1 - ecc.^2);

%make sure it revolves around the sun....
p = [a.*cos(theta); b.*sin(theta)];
p(1,:) = p(1,:) - ecc*a; 
%p = p - p_planet(1:2,1);



u = G*Mass_sun;
v_dir = [-a*sin(theta); b.*cos(theta)]./vecnorm([-a*sin(theta); b.*cos(theta)]);
r_dir = p./vecnorm(p);
vxr_unit = cross([v_dir;zeros(1,N)],[r_dir;zeros(1,N)]);

v = sqrt((1-ecc.^2) * u * a ./ vecnorm(p).^2 ./ vecnorm(vxr_unit).^2) .* v_dir;
for i=1:N
    phi = rand*2*pi;
    rot = [cos(phi), -sin(phi); sin(phi), cos(phi)];
    p(:,i) = rot*p(:,i);
    v(:,i) = rot*v(:,i);
end
p = [p; zeros(1,N)];
v = [v; zeros(1,N)];
%% Output Handling
% -

end

