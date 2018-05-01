clear all; close all;
%Particles in our model;
N = 30;

level_of_awesomeness = 5;
col_index = 1;

%default mass (earth):
defaultMass = 5*10^24;

%Mass of each particle; to be replaced with random; to be replaced by
%actual mass.
Mass = linspace(1,N,N) * defaultMass/10;
Mass(1) = 1000*defaultMass;
Mass = reshape(Mass,[1,N]); %Used for the matrix multiplication in fo

%make a mass combination vector similar as range vector (but then
%multiplied instead of subtracted).
G = 6.67408*10^-11;

%Default range
defaultRange = 5*10^9;

%Make the startposition parmeters (xyz)
X = cos(linspace(1,N,N)/N*2*pi);
Y = sin(linspace(1,N,N)/N*2*pi);
Z = (linspace(1,N,N))*0;

%make the total position vector
p = [X; Y; Z] * defaultRange;
p(:,1) = [0;0;0];

%Defualt speed
defaultSpeed = 5*10^3;

%Make the startspeed parameters;
Vx = cos(linspace(1,N,N)/N*2*pi + pi/2);
Vy = sin(linspace(1,N,N)/N*2*pi + pi/2);
Vz = linspace(1,N,N)*0;

%total speed vector
%v  = [Vx; Vy; Vz] * defaultSpeed;
v = (rand([3 N])-0.5)*2 *defaultSpeed;
v(:,1) = [0;0;0];

% dt = 'stepsize', T = 'total time'
dt = 1000;
T = 1000000000;


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
%define begin energy
E_0 = kin + pot;

%define begin angular momentum
L_0 = AngularMomentum(p,N,Mass,v);

for t = 0:dt:T
    index = index+1;
    
    %later were gonna make some bounds on speed and range, this is needed.
    vOud = v;
    
    %read fo.m first, but keeps track of whether there was a collision.
    c = col(p,Mass,G,N);
    
    %#BUG will crash if multiple collisions in one timestep
    
    %check if the collision vector is empty    
    if max(max(c)) > 0
       m1 = sum(Mass);
       %find indices of collided particles       
       %re-rank the collision indexes
       indices = [mod(find(c),N)'; ceil(find(c)/N)'];
       indices(1,:) = (indices(1,:)==0)*N + indices(1,:);
       indices(:,find(indices(1,:) < indices(2,:))) = [];
       
       
       %ik ga uit van compleet inelastisch.
       %new position is mass centre
       
       M = repmat(Mass,[3 1]);
       p(:,indices(1,:)) = (M(:,indices(1,:)).*p(:,indices(1,:)) + M(:,indices(2,:)).*p(:,indices(2,:)))./(M(:,indices(1,:))+M(:,indices(2,:)));
       
       %use momentum fomulae for new speed
       v(:,indices(1,:)) = (M(:,indices(1,:)).*v(:,indices(1,:)) + M(:,indices(2,:)).*v(:,indices(2,:)))./(M(:,indices(1,:))+M(:,indices(2,:)));
       
       %new mass is um of the masses
       Mass(indices(1,:)) = Mass(indices(1,:)) + Mass(indices(2,:));
       %speed of old particle is 0
       v(:,indices(2,:)) = 0;
       Mass(indices(2,:)) = 0;
       
       %keep track of latest collision index, so that the RMSE error for
       %angular momentum and energy can be better defined.
       col_index = index;
    end
    %calculate the new velocity.
    
    if level_of_awesomeness == 1
        %first order (newton forward)
        a = acc(p,Mass,G,N);
        p = p + v*dt + a * dt^2 / 2;
        v = v + a * dt;
    elseif level_of_awesomeness == 2
        %Runge Kutta 2 conserves angular momentum?
        k1 = dt^2*permute(acc(p + (1/2)*dt*v,Mass,G,N),[3,2,1]);
        p = p + v*dt + k1;
        v = v + k1/dt;
    elseif level_of_awesomeness == 4
        %Runge Kutta 4
        k1 = dt^2*permute(acc(p,Mass,G,N),[3,2,1]);
        k2 = dt^2*permute(acc(p + 0.5*dt*v + 1/8*k1,Mass,G,N),[3,2,1]);
        k3 = dt^2*permute(acc(p + dt*v + .5*k2,Mass,G,N),[3,2,1]);
        p = p + v*dt + 1/6*(k1+2*k2);
        v = v + 1/(6*dt)*(k1+4*k2+k3);
    elseif level_of_awesomeness == 5
        %Runge Kutta 5a see file I (floris) send over whatsapp
        k1 = dt^2*permute(acc(p,Mass,G,N),[3,2,1]);
        k2 = dt^2*permute(acc(p + (1/4)*dt*v + (1/32)*k1,Mass,G,N),[3,2,1]);
        k3 = dt^2*permute(acc(p + (7/10)*dt*v - (7/1000)*k1 + (63/250)*k2,Mass,G,N),[3,2,1]);
        k4 = dt^2*permute(acc(p + dt*v + (2/7)*k1 + (3/14)*k3,Mass,G,N),[3,2,1]);
        p = p + dt*v + (1/14)*k1 + (8/27)*k2 + (25/189)*k4;
        v = v + (1/dt)*((1/14)*k1 + (32/81)*k2 + (250/567)*k3 + (5/54)*k4);
    elseif level_of_awesomeness == 6
        %Runge Kutta 6.
        k1 = dt^2*permute(acc(p,Mass,G,N),[3,2,1]);
        k2 = dt^2*permute(acc(p + 1/4*dt*v + 1/32*k1,Mass,G,N),[3,2,1]);
        k3 = dt^2*permute(acc(p + 1/2*dt*v - k1/24 + k2/6,Mass,G,N),[3,2,1]);
        k4 = dt^2*permute(acc(p + 3/4*dt*v + k1*3/32 + k2/8 + k3/16,Mass,G,N),[3,2,1]);
        k5 = dt^2*permute(acc(p + 3/7*dt*v - k1/14 + k3/7,Mass,G,N),[3,2,1]);
        p = p + dt*v + (7*k1 +24*k2 + 6*k3 + 8*k4)/90;
        v = v + (7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5)/(90*dt);
    end
    
    %als de planeten te ver weg zijn van de oorsprong stop het programma.
    R2 = sqrt(p(1,:).^2 + p(2,:).^2 + p(3,:).^2);
    if min(R2) > defaultRange * 15
        disp('Too far');
        break
    end
    
    [kin,pot] = EnergyTracer(p,N,v,Mass,G);
    
    %make a Total kinetic energy vector for plotting
    T(index) = (kin + pot - E_0) / E_0;
    
    L = AngularMomentum(p,N,Mass,v);
    
    %make a angular momentum vector for plotting
    L_t(index) = (L(3)-L_0(3))/L_0(3);
    
    
    %when plotting too often this can drastically slow down the script. Plotting once every 200 timesteps help speeding this up IFF the plotting is bottlenecking the script
    if mod(index,20) == 0
        subplot(2,2,1) 
        plot(T(max(1,index-5000):end));
        %make the axis nice and kushy
        axis([0 5000 -1 1]);
        %axis([0 5000 min(T(max(1,index-5000):end)) max(T(max(1,index-5000):end))]);
        
        %make the root mean square error of the total energy since the last
        %collision.
        T_RMSE = sqrt(sum((T(col_index:end)-mean(T(col_index(end)))).^2)/index);
        title(T_RMSE);
        drawnow
        
        subplot(2,2,2) 
        plot(L_t(max(1,index-5000):end));
        %make the axis nice and kushy
        axis([[0 5000] [1 1]*round(L_t(end))+[-1 1]]);
        title('Angular momentum(z)')
        drawnow

        subplot(2,2,3); 
        plot(p(1,:),p(2,:),'o');
        
        %Centre around COM;
        CM = COM(Mass,p);
        CM = [CM(1) CM(1) CM(2) CM(2)];
        
        axis(CM + [-1 1 -1 1]*2*defaultRange);
        title(sum(Mass~=0));
        drawnow
    end
    
    %liever niet vol maken want dan heb je nullen
    pos(:,index) = reshape(p,[],1);
end
for i=1:N
    xPos(i,:) = pos(3*(i-1)+1,:);
    yPos(i,:) = pos(3*(i-1)+2,:);
    zPos(i,:) = pos(3*(i-1)+3,:);
end
