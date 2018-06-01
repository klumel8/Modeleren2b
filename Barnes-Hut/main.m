%test results:(100 time steps) (no plotting)
% no Barnes-Hut: (nando's laptop)
% N = 10   -> t = 0.045269
% N = 100  -> t = 0.61684
% N = 1000 -> t = 74.6708
% N = 2000 -> t = 317.9029
% Barnes-hut:(nando's laptop)
% N = 10   -> t = 0.80642
% N = 100  -> t = 0.80225
% N = 1000 -> t = 33.9423
% N = 2000 -> t = 138.5657
powers = 1:4;
N_all = 10.^(powers);
N_all(4) = 2000; %10^4 costs too much time for fast testing

for curr_N = powers
    N = N_all(curr_N);
    
clearvars -except N N_all curr_N test_t; close all;
% clear all; close all;
%Particles in our model;

% N = 100;
G = 6.67408*10^-11; % [Nm^2kg^-2]
defaultRange = 108e9; % [m]

% plotting configuration
fps = 10;

plotting = true;

%integration method
level_of_awesomeness = 7;

%used for plotting
col_index = 1;

%barnes-hut parameters
theta = 1/1.1; 

% Create initial conditions
[Mass, p, v] = initialConditions(defaultRange,N);


% dt = 'stepsize', T = 'total time'

dt = 1000; % in seconds
T = 1e9; % in seconds

%index will later be used to keep track of iterations in order to make a
%plot vector
index = 0;
[kin,pot] = EnergyTracer(p,N,v,Mass,G);
%define begin energy
E_0 = kin + pot;

%define begin angular momentum
L_0 = AngularMomentum(p,N,Mass,v);


%a timer so we dont plot too often and slow down the script
tic;
for t = 0:dt:T
    index = index+1;
    
    %later were gonna make some bounds on speed and range, this is needed.
    vOud = v;
    
    %read fo.m first, but keeps track of whether there was a collision.
    c = col(p,Mass,N);
    
    %#BUG will crash if multiple collisions in one timestep
    
    %check if the collision vector is empty    
    if max(max(c)) > 0
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
    
    if level_of_awesomeness == 1
        %first order (newton forward)
%         a = acc(p,Mass,G,N);
        a = acc_barnes_hut(p,Mass,G,N,theta);
        p = p + v*dt + a * dt^2 / 2;
        v = v + a * dt;
    elseif level_of_awesomeness == 2
        %Runge Kutta 2 conserves angular momentum?
%         k1 = dt^2*acc(p + (1/2)*dt*v,Mass,G,N);
        k1 = dt^2*acc_barnes_hut(p + (1/2)*dt*v,Mass,G,N,theta);
        p = p + v*dt + k1;
        v = v + k1/dt;
    elseif level_of_awesomeness == 4
        %Runge Kutta 4
%         k1 = dt^2*acc(p,Mass,G,N);
%         k2 = dt^2*acc(p + 0.5*dt*v + 1/8*k1,Mass,G,N);
%         k3 = dt^2*acc(p + dt*v + .5*k2,Mass,G,N);
        k1 = dt^2*acc_barnes_hut(p,Mass,G,N,theta);
        k2 = dt^2*acc_barnes_hut(p + 0.5*dt*v + 1/8*k1,Mass,G,N,theta);
        k3 = dt^2*acc_barnes_hut(p + dt*v + .5*k2,Mass,G,N,theta);
        p = p + v*dt + 1/6*(k1+2*k2);
        v = v + 1/(6*dt)*(k1+4*k2+k3);
    elseif level_of_awesomeness == 5
        %Runge Kutta 5a see file I (floris) send over whatsapp
%         k1 = dt^2*acc(p,Mass,G,N);
%         k2 = dt^2*acc(p + (1/4)*dt*v + (1/32)*k1,Mass,G,N);
%         k3 = dt^2*acc(p + (7/10)*dt*v - (7/1000)*k1 + (63/250)*k2,Mass,G,N);
%         k4 = dt^2*acc(p + dt*v + (2/7)*k1 + (3/14)*k3,Mass,G,N);
        k1 = dt^2*acc_barnes_hut(p,Mass,G,N,theta);
        k2 = dt^2*acc_barnes_hut(p + (1/4)*dt*v + (1/32)*k1,Mass,G,N,theta);
        k3 = dt^2*acc_barnes_hut(p + (7/10)*dt*v - (7/1000)*k1 + (63/250)*k2,Mass,G,N,theta);
        k4 = dt^2*acc_barnes_hut(p + dt*v + (2/7)*k1 + (3/14)*k3,Mass,G,N,theta);
        p = p + dt*v + (1/14)*k1 + (8/27)*k2 + (25/189)*k4;
        v = v + (1/dt)*((1/14)*k1 + (32/81)*k2 + (250/567)*k3 + (5/54)*k4);
    elseif level_of_awesomeness == 6
        %Runge Kutta 6.
%         k1 = dt^2*acc(p,Mass,G,N);
%         k2 = dt^2*acc(p + 1/4*dt*v + 1/32*k1,Mass,G,N);
%         k3 = dt^2*acc(p + 1/2*dt*v - k1/24 + k2/6,Mass,G,N);
%         k4 = dt^2*acc(p + 3/4*dt*v + k1*3/32 + k2/8 + k3/16,Mass,G,N);
%         k5 = dt^2*acc(p + 3/7*dt*v - k1/14 + k3/7,Mass,G,N);
        k1 = dt^2*acc_barnes_hut(p,Mass,G,N,theta);
        k2 = dt^2*acc_barnes_hut(p + 1/4*dt*v + 1/32*k1,Mass,G,N,theta);
        k3 = dt^2*acc_barnes_hut(p + 1/2*dt*v - k1/24 + k2/6,Mass,G,N,theta);
        k4 = dt^2*acc_barnes_hut(p + 3/4*dt*v + k1*3/32 + k2/8 + k3/16,Mass,G,N,theta);
        k5 = dt^2*acc_barnes_hut(p + 3/7*dt*v - k1/14 + k3/7,Mass,G,N,theta);
        p = p + dt*v + (7*k1 +24*k2 + 6*k3 + 8*k4)/90;
        v = v + (7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5)/(90*dt);
    elseif level_of_awesomeness == 7
        if t == 0
            %leapfrog initialize acceleration
            %a = acc(p,Mass,G,N);
            a = acc_barnes_hut(p,Mass,G,N,theta);
        end
        %leapfrog
        v = v + a*dt/2;
        p = p + dt*v;
        %a = acc(p,Mass,G,N);
        a = acc_barnes_hut(p,Mass,G,N,theta);
        v = v + a*dt/2;
    end
    
    %fetch the kinetic and potential energy.
    [kin,pot] = EnergyTracer(p,N,v,Mass,G);
    
    %make a Total kinetic energy vector for plotting.
    E_tot(index) = (kin + pot - E_0) / E_0;
    
    %Calculate the Angular momentum
    L = AngularMomentum(p,N,Mass,v);
    
    %make a angular momentum vector for plotting
    L_t(index) = (L(3)-L_0(3))/L_0(3);
    
    %when plotting too often this can drastically slow down the script. Plotting once every 200 timesteps help speeding this up IFF the plotting is bottlenecking the script
    %only plot when 1 == 1, (saves time)

    if mod(t,10*dt) == 0 && plotting%toc > 1/fps && plotting
        figure(1);
        subplot(2,2,1) 
        plot(E_tot);
        axis([max(0,index-5000) index+500 -1 1]);
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', round(xt*dt/31556926,3))
        xlabel('time [years]')
        ylabel('relative magnitude')
        
        %make the root mean square error of the total energy since the last
        %collision.
        E_tot_RMSE = sqrt(sum((E_tot(col_index:end)-mean(E_tot(col_index(end)))).^2)/index);
        title(strcat('RMSE(Energy):',num2str(E_tot_RMSE)));
        
        subplot(2,2,2)       
        plot(L_t);
        title('Angular momentum(z)')
        axis([[max(0,index-5000) index+500] [1 1]*round(L_t(end))+[-1 1]]);
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', round(xt*dt/31556926,3))
        xlabel('time [years]')
        ylabel('relative magnitude')
        

        subplot(2,2,3)
        plot(p(1,2:end),p(2,2:end),'.k','MarkerSize',20); hold on
        plot(p(1,1),p(2,1),'*y', 'MarkerSize',20); hold off
        axis([-1 1 -1 1]*defaultRange);
        axis equal
        title(strcat('N =', " ", num2str(sum(Mass~=0))));
        drawnow
        tic;
    end
end
test_t(curr_N) = toc;
disp('Barnes-Hut')
disp(['N = ', num2str(N)])
disp(['t = ',num2str(test_t(curr_N))])
end

