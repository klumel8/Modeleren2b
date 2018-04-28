function [kin,pot] = EnergyTracer(p,N,v,Mass,G)
%ENERGYTRACER Summary of this function goes here
%   Detailed explanation goes here
    [~,R] = dispVec(p,N);
    pot = -G*(Mass'*Mass)./(R*2);
    pot(pot==-Inf) = 0;
    pot(pot==Inf) = 0;
    pot = nansum(nansum(pot));
    kin = nansum(0.5*Mass.*(v(1,:).^2 + v(2,:).^2 + v(3,:).^2));
    
end

