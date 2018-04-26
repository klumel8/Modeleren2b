clear all; close all;
%Particles in our model;
N = 5;

level_of_awesomeness = 4;

%default mass (earth):
defaultMass = 5*10^24;

%Mass of each particle; to be replaced with random; to be replaced by
%actual mass.
Mass = linspace(1,N,N) * defaultMass;
Mass = reshape(Mass,[1,3]); %Used for the matrix multiplication in fo

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
defaultSpeed = 4*10^2;

%Make the startspeed parameters;
Vx = cos(linspace(1,N,N)/N*2*pi + pi/2);
Vy = sin(linspace(1,N,N)/N*2*pi + pi/2);
Vz = linspace(1,N,N)*0;

%total speed vector
v  = [Vx; Vy; Vz] * defaultSpeed;

% dt = 'stepsize', T = 'total time'
dt = 100;
T = 10000000;


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

for t = 0:dt:T
    index = index+1;
    
    %later were gonna make some bounds on speed and range, this is needed.
    vOud = v;
    
    %read fo.m first, but keeps track of whether there was a collision.
    [~,c] = fo(p,Mass,G,N);
    
    %#BUG will crash if multiple collisions in one timestep
    
    %check if the collision vector is empty
    if max(max(c)) > 0
       
       %find indices of collided particles
       find(c);
       
       %re-rank the collision indexes
       indices = mod(find(c),N);
       indices = (indices==0)*N + indices
       
       
       %ik ga uit van compleet inelastisch.
       %new position is mass centre
       p(:,indices(1)) = (Mass(indices(1))*p(:,indices(1)) + Mass(indices(2))*p(:,indices(2)))/sum(Mass(indices));
       
       %use momentum fomulae for new speed
       v(:,indices(1)) = (Mass(indices(1))*v(:,indices(1)) + Mass(indices(2))*v(:,indices(2)))/sum(Mass(indices));
       
       %new mass is um of the masses
       Mass(indices(1)) = sum(Mass(indices));
       
       %speed of old particle is 0
       v(:,indices(2)) = [0;0;0];
       Mass(indices(2)) = 0;
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
    
    %liever niet vol maken want dan heb je nullen
    pos(:,index) = reshape(p,[],1);
end
for i=1:N
    xPos(i,:) = pos(3*(i-1)+1,:);
    yPos(i,:) = pos(3*(i-1)+2,:);
    zPos(i,:) = pos(3*(i-1)+3,:);
end

hold off
%make the axis stable
minX = min(min(xPos));
maxX = max(max(xPos));
minY = min(min(yPos));
maxY = max(max(yPos));
minZ = min(min(zPos));
maxZ = max(max(zPos));
for t=1:100:size(xPos,2);
    plot3(xPos(:,t),yPos(:,t),zPos(:,t),'*');
    axis([minX maxX minY maxY minZ maxZ]);
    drawnow
end

