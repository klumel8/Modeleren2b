function v = colGas(p, v, Mass, N,dt, MassGas, pGas)
%Calculate which particles, if any, collide. Also include collisions with
%uniformly distributed gass
%input arguments:
%   p       : (3xN) position vector
%   v       : (3xN) velocity vector
%   Mass    : (1xN) Mass vector
%   N       : (scalar) number of particles
%   dt      : (scalar) timestep over which collision is possible
%   MassGas : (scalar) total mass of gas
%   pGas    : (3x ?) Position of gas
%output arguments:
%   Collision: (Nx? logical) the (i,j)th element is 1 if particle i and j
%               collide, 0 if they don't collide. The diagonal is filled
%               with zero's
%required functions(non-standard):
%   dispVec
    rho = 100; % [kgm^-3]
    %gasRho = 0.001; % [kgm^-3]
    %get the distances between the particles
    
    %check if particles collide
    %make the radius of the particle dependent on size
    planetradius = Mass.^(1/3)/(4/3*rho*pi)^(1/3);
    gasRadius = 0.5*(pGas(2,1,1)-pGas(1,1,1))*ones(size(MassGas));
    %do some things to get the linear combination of all added ranges
    radiusCombi = gasRadius + planetradius';
    
    %check whether they collide somewhere
    col_somewhere = (R.^2 - dot(D,DV,3).^2./DV_norm.^2 < radiusCombi.^2);
    %check if the collision is in the interval [0,dt]
    col_time_correct = -dot(D,DV,3)./DV_norm.^2 < dt & -dot(D,DV,3)./DV_norm.^2>0;
    %check if they collide at the begin or endpoint:
    %col_begin_end = R.^2 < planetCombi.^2 | vecnorm(D + DV*dt,2,3).^2 <planetCombi.^2;
    col_begin_end = R.^2 < planetCombi.^2 | sum((D + DV*dt).^2,3) <radiusCombi.^2;
    
    %combine all checks & remove all collisions of a particle with itself
    Collision = ((col_somewhere & col_time_correct) | col_begin_end) & ~eye(N);
    Colsum = sum(Collision,2);
    MassCol = Colsum'*MassGas(1,1);
    v = Mass.*v./(Mass+MassCol);
end
