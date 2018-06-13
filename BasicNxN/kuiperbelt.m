function [p, v, Mass_k] = kuiperbelt(N, p_neptune,trojans)
%kuiperbelt Creates the initial conditions for the simulation
%
%   Syntax:
%
%       [p v]=kuiperbelt(N)
%
%   Input:
%           * N       = "Number of particles to be simulated" (1)[-]
%                        2 = solar system"                    (1)[-]
%           * p_neptune = position of neptune(only used for trojans)
%           * trojans = logical for simulating trojans or not
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

if trojans
    N = 10* N;
end
Mass_k = ones(1,N)*Mass_sun * 1e-10;

% create position and speed vectors
theta = 2*pi*rand(1,N); % create random angles
if trojans
    r = 30.110387*AU;
else
%     r = ((5/2)^(2/3))*(30.110387*AU);
r = (2 + (3.5-2)*rand(1,N))*AU;
% r = (2)^(2/3) * 30.110387*AU;
% r = 39.4*AU;
end

ecc = rand(1,N)*0.01;
% ecc = ones(1,N)*0.1;

%use defualt gonio functions to make a physically valid semi-major/minor
%axis a
a = r;
b = a .* sqrt(1 - ecc.^2);

%make sure it revolves around the sun....
p = [a.*cos(theta); b.*sin(theta)];
p(1,:) = p(1,:) - ecc.*a; 



u = G*Mass_sun;
v_dir = [-a.*sin(theta); b.*cos(theta)]./vecnorm([-a.*sin(theta); b.*cos(theta)]);
r_dir = p./vecnorm(p);
vxr_unit = cross([v_dir;zeros(1,N)],[r_dir;zeros(1,N)]);

v = sqrt((1-ecc.^2) .* u .* a ./ vecnorm(p).^2 ./ vecnorm(vxr_unit).^2) .* v_dir;
for i=1:N
    phi = rand*2*pi;
    if i == 1
        phi = 0;
    end
    rot = [cos(phi), -sin(phi); sin(phi), cos(phi)];
    p(:,i) = rot*p(:,i);
    v(:,i) = rot*v(:,i);
end
p = [p; zeros(1,N)];
v = [v; zeros(1,N)];

if trojans
    wanted_theta = 2*pi/360*[60;-60;120;-120;180];
    difference_theta = 10/360 * 2*pi;
    theta_neptune = atan2(p_neptune(2),p_neptune(1));
    
    trojans_theta_plus = wanted_theta+ difference_theta;
    trojans_theta_min = wanted_theta- difference_theta;
    
    relative_theta = atan2(p(2,:),p(1,:))-theta_neptune;
    relative_theta = relative_theta + (2*pi)* (relative_theta<-pi) - 2*pi*(relative_theta>pi);
    remaining = logical(sum(relative_theta<max(trojans_theta_plus,trojans_theta_min) &...
        relative_theta>min(trojans_theta_plus,trojans_theta_min)));
    
    p = p(:,remaining);
    v = v(:,remaining);
    Mass_k = Mass_k(remaining);
    
    
    
    
    
    
    
end
%% Output Handling
% -

end

