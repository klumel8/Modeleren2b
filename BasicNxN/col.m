function Collision = col(p, v, Mass, N)
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
    rho = 1e3; % [kgm^-3]
    %get the distances between the particles
    [D,R] = dispVec(p,N);    
    [DV,DV_norm] = dispVec(v,N);
        
    %check if particles collide
    %make the radius of the particle dependent on size
    planetradius = Mass.^(1/3)/(4/3*rho*pi)^(1/3);
    
    
    
    
    
    %do some things to get the linear combination of all added ranges
    planetCombi = combvec(planetradius,planetradius);
    planetCombi = planetCombi(1,:) + planetCombi(2,:);
 
    %reshape it back to the standard NxN form
    planetCombi = reshape(planetCombi',[N N]);
    %check whether they collide
    Collision = (R.^2- dot(D,DV,3).^2./DV_norm.^2 < planetCombi).*R>0;
end
