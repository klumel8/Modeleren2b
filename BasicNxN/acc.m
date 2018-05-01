function a = acc(p, Mass, G, N)
    
    p = permute(p,[3,2,1]); %Make the 'xyz' the third dimension. 
    %The first dimension singleton and the 2nd dimension the particles
    repVec = repmat(p,N,1,1); %Repeat the vectors;
    
    D = permute(repVec,[2,1,3])-repVec; %Transpose and subtract the vectors
    %That way subtracting all the different combinations
    
    %Calculate the range between each particle (stored in NxN matrix).
    %R = sqrt(D(:,:,1).^2 + D(:,:,2).^2 + D(:,:,3).^2);
    R = realsqrt(sum(D.^2,3));
    
    %make the mass product;
    massProd = Mass'*Mass; %To make it faster

    
    %Make a simple mathematical expression for everything distance related,
    %later to be used in the force formulae.
    temp= repmat(R, [1,1,3]);
    %gpuTemp = gpuArray(temp);
    %gpuD = gpuArray(D);
    %gpuDisp = gpuD./(gpuTemp.*gpuTemp.*gpuTemp);
    
    dispComp = D./(temp.*temp.*temp);
    
    %calculate the force in each direction (x,y,z).
    %no minus sign because somewhere I misplaced a minus sign (dont worry
    %I anticipated this)
    forceComp = G*repmat(massProd, [1,1,3]).*dispComp;
    
    %calculate the Force of each point on eachother  
    F = nansum(forceComp,1);
    
    %calculate the acceleration in each direction.
    a = F./repmat(Mass,[1,1,3]);
end

