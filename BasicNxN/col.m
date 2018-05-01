function Collision = col(R, Mass)
    
    %check if they collide.
    
    %check if particles collide
    %make the radius of the particle dependent on size
    planetSize = Mass.^(1/3)/(4*pi);
 
    %do some things to get the linear combination of all added ranges
    planetCombi = planetSize+planetSize';
    
    %check whether they collide
    Collision = (R < planetCombi).*R>0;
end
