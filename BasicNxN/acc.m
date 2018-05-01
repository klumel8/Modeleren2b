function a = acc(p, Mass, G, N)
    [D,R] = dispVec(p,N);
    
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

