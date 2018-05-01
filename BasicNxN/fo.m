function [ a, Collision] = fo(p, Mass, G, N)
    [D,R] = dispVec(p,N);
    %make the mass product;
    massProd = Mass'*Mass; %To make it faster
    
    %check if they collide.
    Collision = zeros(N,N);
    
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
    
    %Make a simple mathematical expression for everything distance related,
    %later to be used in the force formulae.
    dispComp = D./(repmat(R, [1,1,3]).^3);
    
    %calculate the force in each direction (x,y,z).
    %no minus sign because somewhere I misplaced a minus sign (dont worry
    %I anticipated this)
    forceComp = G*repmat(massProd, [1,1,3]).*dispComp;
    
    %calculate the Force of each point on eachother
    F = nansum(forceComp,1);
    
    %calculate the acceleration in each direction.
    a = F./repmat(Mass,[1,1,3]);
end
