function [ a, Collision] = fo(p, Mass, G, N)
    %get all combinations of all vectors.
    combiR =  combvec(p,p);
    
    %make the displacement between all vectors.
    D = (combiR(1:3,:) - combiR(4:6,:));
    
    %make the displacement vector in an N x N x 3 matrix. The third
    %dimension stored the x,y,z respectively 1,2,3.
    D = reshape(D',[N N 3]);
    
    %Calculate the range between each particle (stored in NxN matrix).
    R = sqrt(D(:,:,1).^2 + D(:,:,2).^2 + D(:,:,3).^2);
    
    %make the mass product;
    combiM = combvec(Mass,Mass);
    massProd = combiM(1,:) .* combiM(2,:);
    massProd = reshape(massProd',[N N]);
    
    %check if they collide.
    Collision = zeros(N,N);
    colRange = 5*10^3;
    if min(min(R(R~=0))) < colRange
        Collision = (R < colRange).*R>0;
    end
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

