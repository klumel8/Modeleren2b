clear all;
%Particles in our model; compvec
N = 3;

level_of_awesomeness = 4;

%default mass (earth):
defaultMass = 5*10^24;

%Mass of each particle; to be replaced with random; to be replaced by
%actual mass.
Mass = linspace(1,N,N) * defaultMass;
%make a mass combination vector similar as range vector (but then
%multiplied instead of subtracted).
G = 6.67408*10^-11;

%Default range
defaultRange = 10^6;

%Make the startposition parmeters (xyz)
X = cos(linspace(1,N,N)/N*2*pi);
Y = sin(linspace(1,N,N)/N*2*pi);
Z = (linspace(1,N,N));
p = [X; Y; Z] * defaultRange;

%Defualt speed
defaultSpeed = 1*10^4;

%Make the startspeed parameters;
Vx = cos(linspace(1,N,N)/N*2*pi + pi/2);
Vy = sin(linspace(1,N,N)/N*2*pi + pi/2);
Vz = linspace(1,N,N)*0;
v  = [Vx; Vy; Vz] * defaultSpeed;

M = 1/10;
S = 1/10;
T = 10000;
%time step value
step = S*M;

%Total time duration
tVal = T*M;

%pos = zeros(3*N,round(tVal/step));

figure;
hold on;
for i=1:N
    plot(p(1,i),p(2,i),'*');
end

%index will later be used to keep track of iterations in order to make a
%plot vector
index = 0;

for t = 0:step:tVal
    index = index+1;
    
    %later were gonna make some bounds on speed and range, this is needed.
    vOud = v;
    
    [~,c] = fo(p,Mass,G,N);
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
        a = fo(p,Mass,G,N);
        p = p + v*step + a * step^2 / 2;
        v = v + a * step;
    elseif level_of_awesomeness == 2
        %Runge Kutta 2 conserves angular momentum?
        k1 = step^2*squeeze(fo(p + (1/2)*step*v,Mass,G,N))';
        p = p + v*step + k1;
        v = v + k1/step;
    elseif level_of_awesomeness == 4
        %Runge Kutta 4
        k1 = step^2*squeeze(fo(p,Mass,G,N))';
        k2 = step^2*squeeze(fo(p + 0.5*step*v + 1/8*k1,Mass,G,N))';
        k3 = step^2*squeeze(fo(p + step*v + .5*k2,Mass,G,N))';
        p = p + v*step + 1/6*(k1+2*k2);
        v = v + 1/(6*step)*(k1+4*k2+k3);
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
for i=0:N-1
    plot(pos(3*i+1,:),pos(3*i+2,:));
end
