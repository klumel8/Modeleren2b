function [kin,pot] = EnergyTracer(p,N,v,Mass,G)
%ENERGYTRACER calculates the kinetic and the potential energy of the given
%particles
%input arguments:
%   p       : (3xN) position vector
%   N       : (scalar) number of particles
%   v       : (3xN) velocity vector
%   Mass    : (1xN) mass vector
%   G       : (scalar) gravitational constant
%output arguments:
%   kin     : (scalar) total kinetic energy
%   pot     : (scalar) total potential energy
%required functions(non-standard):
%   dispVec

    %first get the distances of each particle to the rest
    [~,R] = dispVec(p,N);
    
    %use the standard formulae for gravitational energy and use r*2 becuase of the symmetry of the R matrix each distance (1 to 2 and 2 to 1) gets counted twice
    pot = -G*(Mass'*Mass)./(R*2);
    
    %filter al the divide by zero potentials
    pot(pot==-Inf) = 0;
    pot(pot==Inf) = 0;
    
    %add them up
    pot = nansum(nansum(pot));
    
    %use the standard formulea for kinetic energy
    kin = nansum(0.5*Mass.*(v(1,:).^2 + v(2,:).^2 + v(3,:).^2));
end
