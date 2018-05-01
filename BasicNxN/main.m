clear; close all;
%Particles in our model; combvec
N = 1000;

level_of_awesomeness = 4; %Which Numerical Integration Method is used

%default mass (earth):
defaultMass = 5*10^24;

%Mass of each particle; to be replaced with random; to be replaced by
%actual mass.
Mass = ones(1,N)* defaultMass;%linspace(1,N,N) * defaultMass;
Mass = reshape(Mass,[1,N]); %Used for the matrix multiplication in fo
%make a mass combination vector similar as range vector (but then
%multiplied instead of subtracted).
G = 6.67408*10^-11;

%Default range
defaultRange = 10^6;

%Make the startposition parmeters (xyz) %Equidistributed over x
X = cos(linspace(1,N,N)/N*2*pi); 
Y = sin(linspace(1,N,N)/N*2*pi); 
Z = (linspace(1,N,N));
p = [X; Y; Z] * defaultRange;

%Default speed
defaultSpeed = 1*10^4;

%Make the startspeed parameters;
Vx = cos(linspace(1,N,N)/N*2*pi + pi/2);
Vy = sin(linspace(1,N,N)/N*2*pi + pi/2);
Vz = linspace(1,N,N)*0;
v  = [Vx; Vy; Vz] * defaultSpeed;


%time dt value
dt = 1;

%Total time duration
T = 100;

%pos = zeros(3*N,round(T/dt));


figure;
hold on;

plot(p(1,:),p(2,:),'*');


%index will later be used to keep track of iterations in order to make a
%plot vector
index = 0;

for t = 0:dt:T
    index = index+1;
    
    %later were gonna make some bounds on speed and range, this is needed.
    vOud = v;
    
    c = col(p,N);
    if max(max(c)) > 0
       find(c);
       indices = mod(find(c),N);
       indices = (indices==0)*N + indices
       %ik ga uit van compleet inelastisch.
       p(:,indices(1)) = (Mass(indices(1))*p(:,indices(1)) + Mass(indices(2))*p(:,indices(2)))/sum(Mass(indices));
       v(:,indices(1)) = (Mass(indices(1))*v(:,indices(1)) + Mass(indices(2))*v(:,indices(2)))/sum(Mass(indices));
       Mass(indices(1)) = sum(Mass(indices));
       v(:,indices(2)) = [0;0;0];
       Mass(indices(2)) = 0;
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
    end
    
    %als de snelheid te snel groter wordt stop dan het programma
    if norm(v(:,:)) > 20*norm(vOud(:,:))
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

for i=0:N-1
    plot3(pos(3*i+1,:),pos(3*i+2,:),pos(3*i+3,:));
end
