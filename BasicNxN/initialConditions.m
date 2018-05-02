function [System] = initialConditions(radius, N)
%initialConditions Creates the initial conditions for the simulation
% Creates a vector with the large mass in the first index, and the other N
% masses spread equally with circular orbits
%   Syntax:
%
%       [System]=initialConditions(radius,N)
%
%   Input:
%           * radius  = "Radius of the orbit of the disc"     (1)[m]
%           * N       = "Number of particles to be simulated" (1)[-]
%   Output:
%           * Mass    = "Mass vector"     (1:(N+1))[kg]
%           * p       = "position vector" (3:(N+1))[m]
%           * v       = "velocity vector" (3:(N+1))[m/s]
%
%   Example(s):
%
%       A system of 1.000.000 particles can be creates by
%               [Mass, p, v] = initialConditions(228e9,1e6)
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
%   Version: v1.0 (02-05-2018)
%

%% Change-log
%
%   v1.0 (02-05-2018)
%       * First operational version.
%

%% Input Handling
%   -

%% Main Function

% Define starting parameters 
Mass_total = (2.6634*10272.66343)e27; % [kg ] total mass in system (no sun)
System(1) = 1.988e30; % [kg] mass of sun 

G = 6.67408*10^-11; % [Nm^2kg^-2]

% calculations
Mass = ones(1,N) * Mass_total / N; % give each partical mass

% create position and speed vectors
theta = 2*pi*rand(1,N); % create random angles
r = radius * ( 0.5 + 0.1*randn(1,N) ); % create normally distributed angles

p = r.*[cos(theta); sin(theta); zeros(1,N)]; % position vector

v_abs = sqrt(G*Mass(1)./r);
v = [sin(v_abs); cos(v_abs); zeros(1,N)];

p(:,1) = [0;0;0]; v(:,1) = [0;0;0]; % pin sun to origin

%% Output Handling
% -

end

