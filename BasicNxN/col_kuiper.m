function Collision = col_kuiper(p,p_k,v,v_k,Mass,dt)
%COL_KUIPER Calculate which kuiperbelt particles, if any, collide with the
%planet(s)
%input arguments:
%   p       : (3xN) position vector of planets
%   p_k     : (3xN_k) position vector of kuiperbelt particles
%   v       : (3xN) velocity vector of planets
%   v_k     : (3xN_k) velocity vector of kuiperbelt particles
%   Mass    : (1xN) Mass vector of planets
%   dt      : (scalar) timestep
%output arguments:
%   Collision: (N_kxN logical) the (i,j)th element is 1 if particle i and j
%               collide, 0 if they don't collide. The diagonal is filled
%               with zero's
%required functions(non-standard):
%   dispVec
    rho = 1e3; % [kgm^-3]
    %get the distances between the particles
%     [D,R] = dispVec(p,N);   
%     [DV,DV_norm] = dispVec(v,N);
    D = permute(p_k,[2,3,1])-permute(p,[3,2,1]);%N_kxNx3
    R = vecnorm(D,2,3);
    DV = permute(v_k,[2,3,1])-permute(v,[3,2,1]);%N_kxNx3
    DV_norm = vecnorm(DV,2,3);
        
    %check if particles collide
    %make the radius of the particle dependent on size
    planetradius = Mass.^(1/3)/(4/3*rho*pi)^(1/3); %(1xN)
    %do some things to get the linear combination of all added ranges
%     planetCombi = planetradius + planetradius';
%     planetCombi = combvec(planetradius,planetradius);
%     planetCombi = planetCombi(1,:) + planetCombi(2,:);
%  
%     %reshape it back to the standard NxN form
%     planetCombi = reshape(planetCombi',[N N]);

    %check whether they collide somewhere
    col_somewhere = (R.^2 - dot(D,DV,3).^2./DV_norm.^2 < planetradius.^2);
    %check if the collision is in the interval [0,dt]
    col_time_correct = -dot(D,DV,3)./DV_norm.^2 < dt & -dot(D,DV,3)./DV_norm.^2>0;
    %check if they collide at the begin or endpoint:
    %col_begin_end = R.^2 < planetCombi.^2 | vecnorm(D + DV*dt,2,3).^2 <planetCombi.^2;
    col_begin_end = R.^2 < planetradius.^2 | sum((D + DV*dt).^2,3) <planetradius.^2;
    
    %combine all checks & remove all collisions of a particle with itself
    Collision = ((col_somewhere & col_time_correct) | col_begin_end);
    
    %old
%     Collision = (R <planetCombi).* R>0;
    
end