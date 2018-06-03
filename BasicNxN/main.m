clear; close all;
%Particles in our model;
% 1 = early solar system
% 2 = solar system and Kuyper belt
% 3 = sphere
type = 2;

%integration method
%1: newton forward
%2,4-6: runge kutta 
%7: leapfrog
int_met = 7;
%use barnes hut
barnes_hut = false;

% universal parameters
G = 6.67408*10^-11; % [Nm^2kg^-2]
AU = 1.49597871e11; % [m]

if type == 1 % early solar system
    defaultRange = 5*AU; % [m]
    N = 1e2;
    dt = 3600*24; % in seconds (dt = 1 day)
    T = 5e10; % in seconds
    [Mass, p, v, N] = initialConditions(defaultRange,N,1);
    
end

if type == 2 % solar system and Kuyper belt
    defaultRange = 5e12; % [m]
    N = 1; % Dummy variable
    N_k = 1e3; % particles in kuiper belt
    dt = 3600*24*7*52; % in seconds 
    T = 1e12; % in seconds
    [Mass, p, v, N] = initialConditions(defaultRange,N,2);
    [p_k, v_k] = kuiperbelt(N_k);
end

% plotting configuration
fps = 10;
plot_system = true;     %plot the particle system
plot_ecc_a = true;      %plot eccentricity vs semi major axis
plot_ang_mom = true;    %plot the angular momentum
plot_momentum = false;   %plot the momentum, relative to jupiter(only for type ==2)
plotting = true;        %plot anything at all

%remove the particles every [remove_index] timesteps
remove_index = 10;
removing = true;

% colision index used for plotting
colision_index = 1;

%index will later be used to keep track of iterations in order to make a
%plot vector
index = 0;
[kin,pot] = EnergyTracer(p,N,v,Mass,G);
%define begin energy
E_0 = kin + pot;

%define begin angular momentum
L_0 = AngularMomentum(p,N,Mass,v);

% a timer so we dont plot too often and slow down the script
tic;
for t = 0:dt:T
    index = index+1;
    
    %remove particles every [remove_index] timesteps:
    if mod(index,remove_index) == 0 && removing 
        %to be removed from: p,v,Mass,N
        %only select the indices which wont be removed:
        staying_indices = find(Mass ~= 0 & ~isnan(Mass));
        p = p(:,staying_indices); 
        v = v(:,staying_indices);
        Mass = Mass(staying_indices);
        N = numel(staying_indices);
        if int_met == 7
            a = a(:,staying_indices);
        end
    end
    
    %later were gonna make some bounds on speed and range, this is needed.
    vOud = v;
    
    %remove particles which are too far with a too high speed:
    remove_crit = vecnorm(p)>10*defaultRange;
    if any(remove_crit)
        indices = find(remove_crit);
        p(:,indices) = p(:,indices)./20;%bring particle back in the system
        v(:,indices) = repmat([0;0;0],1,numel(indices));%set velocity to 0
        Mass(indices) = 0;%set mass to 0
    end
    
    %read fo.m first, but keeps track of whether there was a collision.
    c = col(p,v,Mass,N,dt);
    
    %#BUG will crash if multiple collisions in one timestep
    %check if the collision vector is empty    
    if max(max(c)) > 0
       %find indices of collided particles       
       %re-rank the collision indexes
       indices = [mod(find(c),N)'; ceil(find(c)/N)'];
       indices(1,:) = (indices(1,:)==0)*N + indices(1,:);
       indices(:,find(indices(1,:) > indices(2,:))) = [];
       %this way, if a particle collides with the sun, the sun is still the first particle
        
       
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
       colision_index = index;
    end
    
    if int_met == 1
        %first order (newton forward)
        a = acc(p,Mass,G,N);
        p = p + v*dt + a * dt^2 / 2;
        v = v + a * dt;
    elseif int_met == 2
        %Runge Kutta 2 conserves angular momentum?
        k1 = dt^2*acc(p + (1/2)*dt*v,Mass,G,N);
        p = p + v*dt + k1;
        v = v + k1/dt;
    elseif int_met == 4
        %Runge Kutta 4
        k1 = dt^2*acc(p,Mass,G,N);
        k2 = dt^2*acc(p + 0.5*dt*v + 1/8*k1,Mass,G,N);
        k3 = dt^2*acc(p + dt*v + .5*k2,Mass,G,N);
        p = p + v*dt + 1/6*(k1+2*k2);
        v = v + 1/(6*dt)*(k1+4*k2+k3);
        
    elseif int_met == 5
        %Runge Kutta 5a see file I (floris) send over whatsapp
        k1 = dt^2*acc(p,Mass,G,N);
        k2 = dt^2*acc(p + (1/4)*dt*v + (1/32)*k1,Mass,G,N);
        k3 = dt^2*acc(p + (7/10)*dt*v - (7/1000)*k1 + (63/250)*k2,Mass,G,N);
        k4 = dt^2*acc(p + dt*v + (2/7)*k1 + (3/14)*k3,Mass,G,N);
        p = p + dt*v + (1/14)*k1 + (8/27)*k2 + (25/189)*k4;
        v = v + (1/dt)*((1/14)*k1 + (32/81)*k2 + (250/567)*k3 + (5/54)*k4);
    elseif int_met == 6
        %Runge Kutta 6.
        k1 = dt^2*acc(p,Mass,G,N);
        k2 = dt^2*acc(p + 1/4*dt*v + 1/32*k1,Mass,G,N);
        k3 = dt^2*acc(p + 1/2*dt*v - k1/24 + k2/6,Mass,G,N);
        k4 = dt^2*acc(p + 3/4*dt*v + k1*3/32 + k2/8 + k3/16,Mass,G,N);
        k5 = dt^2*acc(p + 3/7*dt*v - k1/14 + k3/7,Mass,G,N);
        p = p + dt*v + (7*k1 +24*k2 + 6*k3 + 8*k4)/90;
        v = v + (7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5)/(90*dt);
    elseif int_met == 7
        if t == 0
            %initialize acceleration for leapfrog
            a = acc(p,Mass,G,N);
        end
        %leapfrog
        v = v + dt/2*a;
        p = p + dt*v;
        a = acc(p,Mass,G,N);
        v = v + a*dt/2;
        
        if type == 2
            if t == 0
                %initialize acceleration for leapfrog
                a_k = kuiperacc(p,p_k,Mass);
            end
            %leapfrog
            v_k = v_k + dt/2*a_k;
            p_k = p_k + dt*v_k;
            a_k = kuiperacc(p,p_k,Mass);
            v_k = v_k + a_k*dt/2;            
        end
    end
    
    %fetch the kinetic and potential energy.
    [kin,pot] = EnergyTracer(p,N,v,Mass,G);
    
    %make a Total kinetic energy vector for plotting.
    E_tot(index) = (kin + pot - E_0) / E_0;
    
    %Calculate the Angular momentum
    L = AngularMomentum(p,N,Mass,v);
    
    %make a angular momentum vector for plotting
    L_t(index) = (L(3)-L_0(3))/L_0(3);
    
    if type ~= 2
        %make a momentum vector for plotting (only the norm)
        momentum(:,:,index) = Mass.*v; % momentum of all particles (3xNxtime)
        momentum_norm = vecnorm(nansum(momentum,2),2,1); %(1x1xtime)
        rel_momentum = momentum_norm./vecnorm(momentum(:,end-3,1),2,1); %momentum relative to jupiter
        rel_momentum = permute(rel_momentum,[3,2,1]);
    end
    
    [ecc, semi_m_axis] = eccentricity_sma(p,v,Mass);
    ecc = vecnorm(ecc,2,1)';
    semi_m_axis = semi_m_axis';

    %when plotting too often this can drastically slow down the script. Plotting once every 200 timesteps help speeding this up IFF the plotting is bottlenecking the script
    %only plot when 1 == 1, (saves time)
    if toc > 1/fps && plotting
        figure(1)
        if plot_ecc_a
            %eccentricity vs semi-major axis: 
            subplot(2,2,1) 
            plot(semi_m_axis(2:end),ecc(2:end),'.')
            axis([0,1.1*defaultRange,0, 0.5])
            title(['time: ',num2str(round(t/31556926,1)),' y'])
            ylabel('$\varepsilon$','Interpreter','Latex')
            xlabel('a[m]')
        end
        if plot_ang_mom
            %angular momentum
            subplot(2,2,2)       
            plot(L_t);
            title('Angular momentum(z)')
            axis([[max(0,index-5000) index+500] [1 1]*round(L_t(end))+[-1 1]]);
            xt = get(gca, 'XTick');
            set(gca, 'XTick', xt, 'XTickLabel', round(xt*dt/31556926,1))
            xlabel('time [years]')
            ylabel('relative magnitude')
        end
        if plot_system
            T_neptune = 60182*3600*24; % seconds
            omega_neptune = 2*pi/T_neptune;
            A = [cos(omega_neptune*t), sin(omega_neptune*t);...
                -sin(omega_neptune*t), cos(omega_neptune*t) ];
            plot_p(1:2,:) = A*p(1:2,:);
            plot_p_k(1:2,:) = A*p_k(1:2,:);
            %particle system
            subplot(2,2,3)
            plot(plot_p(1,2:end),plot_p(2,2:end),'.k','MarkerSize',20); hold on
            plot(plot_p(1,1),plot_p(2,1),'*y', 'MarkerSize',20); hold on
            axis([-1 1 -1 1]*defaultRange*1.1);
            title(strcat('N =', " ", num2str(sum(Mass~=0)-1)));
            
            %plot kuiperbelt if type == 2
            if type == 2
                plot(plot_p_k(1,:),plot_p_k(2,:),'.r','MarkerSize',2); hold on
                axis([-1 1 -1 1]*50*AU);
            end
            hold off;
        end

        if type == 2
            if plot_momentum
                subplot(2,2,4)
%                 plot(rel_momentum);
                title('Rel momentum(norm), rel to jupiter')
                axis([max(0,index-5000) index+500 -0.1 1]);

                xt = get(gca, 'XTick');
                set(gca, 'XTick', xt, 'XTickLabel', round(xt*dt/31556926,1))
                xlabel('time [years]')
                ylabel('relative magnitude')
            end
        end
        drawnow
        tic;
    end
end
