function Collision = col(p, v, Mass, N,dt)
%Calculate which particles, if any, collide
%input arguments:
%   p       : (3xN) position vector
%   v       : (3xN) velocity vector
%   Mass    : (1xN) Mass vector
%   N       : (scalar) number of particles
%output arguments:
%   Collision: (NxN logical) the (i,j)th element is 1 if particle i and j
%               collide, 0 if they don't collide. The diagonal is filled
%               with zero's
%required functions(non-standard):
%   dispVec
    rho = 100; % [kgm^-3]
    %get the distances between the particles
    [D,R] = dispVec(p,N);   
    [DV,DV_norm] = dispVec(v,N);
        
    %check if particles collide
    %make the radius of the particle dependent on size
    planetradius = Mass.^(1/3)/(4/3*rho*pi)^(1/3);
    if type == 1 || type == 2 || type == 5
        planetradius(1) = 6.9e8; %radius of the sun
    end
    
    %do some things to get the linear combination of all added ranges
    planetCombi = planetradius + planetradius';
    
    %check whether they collide somewhere
    col_somewhere = (R.^2 - dot(D,DV,3).^2./DV_norm.^2 < planetCombi.^2);
    %check if the collision is in the interval [0,dt]
    col_time_correct = -dot(D,DV,3)./DV_norm.^2 < dt & -dot(D,DV,3)./DV_norm.^2>0;
    %check if they collide at the begin or endpoint:
    %col_begin_end = R.^2 < planetCombi.^2 | vecnorm(D + DV*dt,2,3).^2 <planetCombi.^2;
    col_begin_end = R.^2 < planetCombi.^2 | sum((D + DV*dt).^2,3) <planetCombi.^2;
    
    %combine all checks & remove all collisions of a particle with itself
    Collision = ((col_somewhere & col_time_correct) | col_begin_end) & ~eye(N);
    
end
