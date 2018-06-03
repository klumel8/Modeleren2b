function [p, v] = kuiperbelt(N)
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
r_low = 42*AU;
r_high = 48*AU;

% create position and speed vectors
theta = 2*pi*rand(1,N); % create random angles
% r = r_low + (r_high-r_low).*rand(1,N); % create uniformly distributed radii
r = (3/2)^(2/3)*4495e9;
p = r.*[cos(theta); sin(theta); zeros(1,N)]; % position vector

v_abs = sqrt(G*Mass_sun ./ r);
v = v_abs .* [-sin(theta); cos(theta); zeros(1,N)];

%% Output Handling
% -

end
