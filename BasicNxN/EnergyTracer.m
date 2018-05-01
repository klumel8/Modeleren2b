function [kin,pot] = EnergyTracer(R,v,Mass,G)
%ENERGYTRACER Summary of this function goes here
%   The goal of this function is to keep track fo the total energy of the system

    
    %use the standard formulae for gravitational energy and use r*2 becuase of the symmetry of the R matrix each distance (1 to 2 and 2 to 1) gets counted twice
    potMat = -G*(Mass'*Mass)./(R*2);
    
    %filter al the divide by zero potentials
    potMat(potMat==-Inf) = 0;
    potMat(potMat==Inf) = 0;
    
    %add them up
    pot = nansum(nansum(potMat));
    
    %use the standard formulea for kinetic energy
    kin = nansum(0.5*Mass.*sum(v.^2,1)); 
end
