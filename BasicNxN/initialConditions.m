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

if type == 1 %random particles
    % calculations
    Mass = ones(1,N) * Mass_total / (N-1); % give each partical mass
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
   
if type == 2 %kuiperbelt
% %          Sun Mercury Venus Earth Mars  Jupiter Saturn Uranus Neptune
%     Mass= [0.1 0.330   4.87  5.97  0.64  1898    568    68.6   102   ] * 10^24;   % kg
%     r   = [0.1 57.9    108.2 149.6 227.9 778.6   1433.5 2872.5 4495.1] * 10^9;    % m

%          Sun Jupiter Saturn Uranus Neptune
%     Mass= [0.1 1898    568    68.6   102   ] * 10^24;   % kg
%     r   = [0.1 778.6   1433.5 2872.5 4495.1] * 10^9;    % m
    
    %for richardson error estimation
%     Mass = [0.1, 1898] * 10^24;
%     r   = [0.1 778.6 ]* 10^9;


%          Sun Neptune
    Mass= [0.1 102] * 10^24;   % kg % not anymore:mass of neptune is 4* 100 times higher than normal
    r   = [0.1 4495.1] * 10^9;    % m
    N   = length(Mass);
    
    Mass(1) = Mass_sun;
    %for richardson error estimation:
%     theta = 2*pi*linspace(0,1,numel(r));

    theta = 2*pi*[0 0.25]; % create random angles
      
    % neptune information
    ecc = [0,0.009456];
    a   = r;
    b = a.*sqrt(1-ecc.^2);
    r = (a.*(1 - ecc.^2) ) ./ (1 + ecc.*cos(theta));
    v_abs = sqrt(G * Mass(1) * (2./r - 1./a));
    p = [a.*cos(theta); b.*sin(theta); zeros(1,N) ];
    v = v_abs .* [a.*-sin(theta); b.*cos(theta); zeros(1,N) ]./(vecnorm([a.*-sin(theta); b.*cos(theta); zeros(1,N) ]));
    
    % pin sun to origin
    p(:,1) = [0;0;0]; v(:,1) = [0;0;0]; 
    momentum = nansum(Mass .* v,2);
    v(:,1) = -momentum / Mass(1); %give sun position and velocity to make velocitiy of CoM 0
    p(:,1) = -nansum(Mass.*p,2)/Mass(1);
end

if type == 3 %sphere
    Mass_total = 2.6634e32;
    Mass = rand(1,N)/N; % give each partical mass
    Mass = Mass*Mass_total/(sum(Mass));
    % create position and speed vectors
    theta   = 2*pi*rand(1,N); % create random angles
    phi     =   pi*rand(1,N);
    r = radius * ( rand(1,N) ); % create normally distributed angles

    p = r.*[cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)]; % position vector

    
    %v = 10* randn(3,N) + sqrt(G*Mass./(100*radius)).*[-sin(theta);cos(theta);zeros(1,N)];
    v = 10*randn(3,N) + 10^2*[-sin(theta);cos(theta);zeros(1,N)];
%     v = zeros(3,N);
%     v_abs = sqrt(G*Mass(1)./(r.*abs(sin(phi))));
%     v = v_abs .* [-sin(theta); cos(theta); zeros(1,N)];
    

end
if type  == 4 %test to see whether calculations are correct
   Mass = [2,1]*Mass_total/3;
   p = [-0.5,0.5;0,0;0,0]*radius;
   v = [0,0;-1/2,1;0,0]*sqrt(G*Mass(2)./(2*radius));
   
end
if type == 5 %normal solar system
    %          Sun Mercury Venus Earth Mars  Jupiter Saturn Uranus Neptune
    Mass= [0.1 0.330   4.87  5.97  0.64  1898    568    68.6   102   ] * 10^24;   % kg
    r   = [0.1 57.9    108.2 149.6 227.9 778.6   1433.5 2872.5 4495.1] * 10^9;    % m

    N   = length(Mass);
    
    Mass(1) = Mass_sun;
    theta = 2*pi*rand(1,N); % create random angles
    p = r.*[cos(theta); sin(theta); zeros(1,N)]; % position vector

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

