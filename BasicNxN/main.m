clear all; close all;
%Particles in our model;
N = 5;

level_of_awesomeness = 4;

%default mass (earth):
defaultMass = 5*10^24;

%Mass of each particle; to be replaced with random; to be replaced by
%actual mass.
Mass = linspace(1,N,N) * defaultMass;
Mass = reshape(Mass,[1,N]); %Used for the matrix multiplication in fo

%make a mass combination vector similar as range vector (but then
%multiplied instead of subtracted).
G = 6.67408*10^-11;

%Default range
defaultRange = 10^9;

%Make the startposition parmeters (xyz)
X = cos(linspace(1,N,N)/N*2*pi);
Y = sin(linspace(1,N,N)/N*2*pi);
Z = (linspace(1,N,N));

%make the total position vector
p = [X; Y; Z] * defaultRange;

%Defualt speed
defaultSpeed = 1*10^2;

%Make the startspeed parameters;
Vx = cos(linspace(1,N,N)/N*2*pi + pi/2);
Vy = sin(linspace(1,N,N)/N*2*pi + pi/2);
Vz = linspace(1,N,N)*0;

%total speed vector
v  = [Vx; Vy; Vz] * defaultSpeed;

% dt = 'stepsize', T = 'total time'
dt = 1000;
T = 100000000;


%pos = zeros(3*N,round(T/dt));

%some old plot code
% figure;
% hold on;
% for i=1:N
%     plot(p(1,i),p(2,i),'*');
% end

%index will later be used to keep track of iterations in order to make a
%plot vector
index = 0;
[kin,pot] = EnergyTracer(p,N,v,Mass,G);
E_0 = kin + pot;
for t = 0:dt:T
    index = index+1;
    
    %later were gonna make some bounds on speed and range, this is needed.
    vOud = v;
    
    %read fo.m first, but keeps track of whether there was a collision.
    [~,c] = fo(p,Mass,G,N);
    
    %#BUG will crash if multiple collisions in one timestep
    
    %check if the collision vector is empty    
    if max(max(c)) > 0
       m1 = sum(Mass);
       %find indices of collided particles       
       %re-rank the collision indexes
       %indices are the [a,b] coordinates of the NxN matrix which describes the particles which collise, a equals the column number and b the row number.
       indices = [mod(find(c),N)'; ceil(find(c)/N)'];
       
       %mod(5,5) returns 0, but this should be 5, so we fix this with this easy step.
       indices(1,:) = (indices(1,:)==0)*N + indices(1,:);
       
       %here we select only the elements above the diagonal becuase else we would have duplicates.
       indices(:,find(indices(1,:) < indices(2,:))) = [];
       
       
       %ik ga uit van compleet inelastisch.
       %new position is mass centre
       
       %the mass centre is the new centre
       M = repmat(Mass,[3 1]);
       p(:,indices(1,:)) = (M(:,indices(1,:)).*p(:,indices(1,:)) + M(:,indices(2,:)).*p(:,indices(2,:)))./(M(:,indices(1,:))+M(:,indices(2,:)));
       
       %use momentum fomulae for new speed
       v(:,indices(1,:)) = (M(:,indices(1,:)).*v(:,indices(1,:)) + M(:,indices(2,:)).*v(:,indices(2,:)))./(M(:,indices(1,:))+M(:,indices(2,:)));
       
       %new mass is um of the masses
       Mass(indices(1,:)) = Mass(indices(1,:)) + Mass(indices(2,:));
       %speed of old particle is 0
       v(:,indices(2,:)) = 0;
       Mass(indices(2,:)) = 0;
       
       %show a message if there was significant (non rounding error) Mass loss.
       if abs(m1 - sum(Mass)) > 10^13
           m1 - sum(Mass)
           indices
       end
    end
    %calculate the new velocity.
    
    if level_of_awesomeness == 1
        %first order (newton forward)
        a = fo(p,Mass,G,N);
        p = p + v*dt + a * dt^2 / 2;
        v = v + a * dt;
    elseif level_of_awesomeness == 2
        %Runge Kutta 2 conserves angular momentum?
        k1 = dt^2*permute(fo(p + (1/2)*dt*v,Mass,G,N),[3,2,1]);
        p = p + v*dt + k1;
        v = v + k1/dt;
    elseif level_of_awesomeness == 4
        %Runge Kutta 4
        k1 = dt^2*permute(fo(p,Mass,G,N),[3,2,1]);
        k2 = dt^2*permute(fo(p + 0.5*dt*v + 1/8*k1,Mass,G,N),[3,2,1]);
        k3 = dt^2*permute(fo(p + dt*v + .5*k2,Mass,G,N),[3,2,1]);
        p = p + v*dt + 1/6*(k1+2*k2);
        v = v + 1/(6*dt)*(k1+4*k2+k3);
    end
    
    %als de snelheid te snel groter wordt stop dan het programma
    if norm(v(:,:)) > 20*norm(vOud(:,:));
        disp('Too fast')
        break
    end
    
    %als de planeten te ver weg zijn van de oorsprong stop het programma.
    R2 = sqrt(p(1,:).^2 + p(2,:).^2 + p(3,:).^2);
    if min(R2) > defaultRange * 5
        disp('Too far');
        break
    end
    
    %Make a Total energy vector for plotting
    [kin,pot] = EnergyTracer(p,N,v,Mass,G);
    T(index) = kin + pot - E_0;
    
    %when plotting too often this can drastically slow down the script. Plotting once every 200 timesteps help speeding this up IFF the plotting is bottlenecking the script
    if mod(index,200) == 0
        subplot(1,2,1) 
        plot(T(max(1,index-5000):end));
        %make the axis nice and kushy
        axis([0 5000 min(T(max(1,index-5000):end)) max(T(max(1,index-5000):end))]);
        title('total energy')
        drawnow

        subplot(1,2,2) 
        plot(p(1,:),p(2,:),'o');
        axis([-1 1 -1 1]*10*defaultRange);
        title('Movement')
        drawnow
    end
end
