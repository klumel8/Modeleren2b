function [Mass, p, v, N] = initialConditions(radius, N, type)
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
%           * type    = "1 = loose particles only,
%                        2 = solar system"                    (1)[-]
%   Output:
%           * Mass    = "Mass vector"     (1:N)[kg]
%           * p       = "position vector" (3:N)[m]
%           * v       = "velocity vector" (3:N)[m/s]
%           * N       = "Number of particles to be simulated" (1)[-]
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
Mass_total = 2.6634e27; % [kg ] total mass in system (no sun)
Mass_sun = 1.988e30; % [kg] mass of sun 
G = 6.67408*10^-11; % [Nm^2kg^-2]
AU = 1.49597871e11;% [m]
mu = G*Mass_sun; % standard gravitational parameter

if type == 1
    % calculations
    Mass = ones(1,N) * Mass_total / N; % give each partical mass
    Mass(1) = Mass_sun;

    % create position and speed vectors
    theta = 2*pi*rand(1,N); % create random angles
    AU = 1.49597871e11;% [m]
    r_low = 1*AU;
    r_high = 5*AU;
    r = r_low + (r_high - r_low)*rand(1,N); % create normally distributed angles

    p = r.*[cos(theta); sin(theta); zeros(1,N)]; % position vector

    v_abs = sqrt(G*Mass(1)./r);
    v = v_abs .* [-sin(theta); cos(theta); zeros(1,N)];

    p(:,1) = [0;0;0]; v(:,1) = [0;0;0]; % pin sun to origin
    momentum = nansum(Mass .* v,2);
    v(:,1) = -momentum / Mass(1); %give sun position and velocity to make velocitiy of CoM 0
    p(:,1) = -nansum(Mass.*p,2)/Mass(1);
end
   
if type == 2
    % Create planets
    %        Sun Neptune   Pluto                   Jupiter
%     Mass  = [0   102       0.013                   1898 ] * 10^24;   % kg
%     r     = [0.1 30.110387 ((3/2)^(2/3))*30.110387 5.2  ] * AU;      % AU
%     ecc   = [0   0.009456  0.25                    0.05 ];
%     theta = [0   0         180/360                 0    ] * 2*pi;    % radians
    
    Mass  = [0   102       0.013                   ] * 10^24;   % kg
    r     = [0.1 30.110387 ((3/2)^(2/3))*30.110387 ] * AU;      % AU
    ecc   = [0   0.009456  0.25                    ];
    theta = [0   0         180/360                 ] * 2*pi;    % radians
    
    Mass(1) =  Mass_sun;
    N     = length(Mass);
    
    % Define semi major axis and semi minor axis
    a = r;
    b = a .* sqrt(1 - ecc.^2);

    % Create position vector and place focus in origin
    p = [a.*cos(theta); b.*sin(theta)];
    p(1,:) = p(1,:) - ecc.*a; 

    % Calculate speed and distance to origin
    v_abs = vecnorm([-a.*sin(theta); b.*cos(theta)]);
    p_abs = vecnorm(p);
    
    % Create direction vectors position and velocity
    v_dir = [-a.*sin(theta); b.*cos(theta)]./v_abs;
    r_dir = p./p_abs;
    
    % Cross product v and r to calculate total velocity
    vxr_unit = cross([v_dir;zeros(1,N)],[r_dir;zeros(1,N)]);
    v = sqrt((1-ecc.^2) .* mu .* a ./ vecnorm(p).^2 ./ vecnorm(vxr_unit).^2) .* v_dir;
    
    % create z direction
    p = [p; zeros(1,N)];
    v = [v; zeros(1,N)];
    
    % pin sun to origin
    p(:,1) = [0;0;0]; v(:,1) = [0;0;0]; 
    momentum = nansum(Mass .* v,2);
    v(:,1) = -momentum / Mass(1); %give sun position and velocity to make velocitiy of CoM 0
    p(:,1) = -nansum(Mass.*p,2)/Mass(1);
end

if type == 3
    Mass = ones(1,N) * Mass_total / N; % give each partical mass

    % create position and speed vectors
    theta   = 2*pi*rand(1,N); % create random angles
    phi     =   pi*rand(1,N);
    r = radius * ( 0.5 + 0.1*randn(1,N) ); % create normally distributed angles

    p = r.*[cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)]; % position vector

    v_abs = sqrt(G*Mass(1)./r);
    v = v_abs .* [-sin(theta); cos(theta); zeros(1,N)];

    p(:,1) = [0;0;0]; v(:,1) = [0;0;0]; % pin sun to origin   
    momentum = nansum(Mass .* v,2);
    v(:,1) = -momentum / Mass(1);
    p(:,1) = -nansum(Mass.*p,2)/Mass(1);
end
%% Output Handling
% -

end

