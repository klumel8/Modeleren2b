function Collision = col(p, N)
    [~,R] = dispVec(p,N);
    
    %check if they collide.
    
        
    %check if particles collide
    %make the radius of the particle dependent on size
    planetSize = Mass.^(1/3)/(4*pi);
 
    %do some things to get the linear combination of all added ranges
    planetCombi = combvec(planetSize,planetSize);
    planetCombi = planetCombi(1,:) + planetCombi(2,:);
 
    %reshape it back to the standard NxN form
    planetCombi = reshape(planetCombi',[N N]);
    %check whether they collide
    Collision = (R < planetCombi).*R>0;
end
