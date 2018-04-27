function [kin,pot] = EnergyTracer(p,N,v,Mass,G)
%ENERGYTRACER Summary of this function goes here
%   The goal of this function is to keep track fo the total energy of the system

    %first get the distances of each particle to the rest
    [~,R] = dispVec(p,N);
    
    %use the standard formulae for gravitational energy and use r*2 becuase of the symmetry of the R matrix each distance (1 to 2 and 2 to 1) gets counted twice
    pot = -G*(Mass'*Mass)./(R*2);
    
    %filter al the divide by zero potentials
    pot(pot==-Inf) = 0;
    pot(pot==Inf) = 0;
    
    %add them up
    pot = nansum(nansum(pot))
    
    %use the standard formulea for kinetic energy
    kin = nansum(0.5*Mass.*(v(1,:).^2 + v(2,:).^2 + v(3,:).^2)); 
end
