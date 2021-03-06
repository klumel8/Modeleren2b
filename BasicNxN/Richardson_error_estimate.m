%different main for calculating richardson error estimation(REE)
%when running this, change the theta to a linspace(commented) and the
%planets to 2 in initialConditions

clear all; close all;
%all integration methods where REE is calculated for:
%integration method
%1: newton forward
%2,4-6: runge kutta (paper)
%7: leapfrog
%8: runge kutta 4, first order
int_methods = [4,7];%[1,2,4,5,6,7,8];

for int_met = int_methods
    start_dt = 3600*24*7*52*50;
    T = start_dt;
    saved_var = [];
    P = 13;
    xCells = cell(P,1);
    for power = 1:P
        dt = start_dt*2^(-(power));
        %Particles in our model;
        % 1 = early solar system
        % 2 = solar system and Kuyper belt
        % 3 = sphere
        type = 2;
        xPower = zeros(2,T/dt+1);
        %use barnes hut
        barnes_hut = false;
        
        % universal parameters
        G = 6.67408*10^-11; % [Nm^2kg^-2]
        AU = 1.49597871e11; % [m]
        rng(121)
        if type == 1 % early solar system
            defaultRange = 5*AU; % [m]
            N = 2;
            %     dt = 3600*24; % in seconds (dt = 1 day)
            %     T = 5e10; % in seconds
            [Mass, p, v, N] = initialConditions(defaultRange,N,1);
        end
        
        
        if type == 2 % solar system and Kuyper belt
            defaultRange = 5e12; % [m]
            N = 1; % Dummy variable
            N_k = 0; % particles in kuiper belt
            %     dt = 26*3600*24*7; % in seconds (dt = 2 weeks)
            %     T = 1e12; % in seconds
            [Mass, p, v, N] = initialConditions(defaultRange,N,2);
            [p_k, v_k] = kuiperbelt(N_k);
        end
        p(2,1)
        
        begin_p = p;
        % plotting configuration
        fps = 10;
        plot_system = true;     %plot the particle system
        plot_ecc_a = true;      %plot eccentricity vs semi major axis
        plot_ang_mom = true;    %plot the angular momentum
        plot_momentum = true;   %plot the momentum, relative to jupiter(only for type ==2)
        plotting = false;        %plot anything at all
        
        %remove the particles every [remove_index] timesteps
        remove_index = 10;
        removing = false;
        
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
            xPower(:,t/dt+1) = [p(1,2);t];
            index = index+1;
            
            %remove particles every [remove_index] timesteps:
            if mod(index,remove_index) == 0 && removing
                %to be removed from: p,v,Mass,N
                %only select the indices which wont be removed:
                staying_indices = find(Mass ~= 0 & ~isnan(Mass));
                disp('test')
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
                %         kv1 = acc(p,Mass,G,N);
                %         kr1 = v;
                %         kv2 = acc(p+1/2*dt*kr1,Mass,G,N);
                %         kr2 = v.*kv1*dt/2;
                %         kv3 = acc(p+1/2*dt*kr2,Mass,G,N);
                %         kr3 = v.*kv2*dt/2;
                %         kv4 = acc(p+kr3*dt,Mass,G,N);
                %         kr4 = v.*kv3*dt;
                %         v = v + dt/6*(kv1+2*kv2+2*kv3+kv4);
                %         p = p + dt/6*(kr1+2*kr2+2*kr3+kr4);
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
                    %             if t == 0
                    %                 %initialize acceleration for leapfrog
                    %                 a = acc(p_k,Mass,G,N);
                    %             end
                    %             %leapfrog
                    %             v = v + dt/2*a;
                    %             p = p + dt*v;
                    %             a = acc(p,Mass,G,N);
                    %             v = v + a*dt/2;
                end
            elseif int_met == 8
                kv1 = dt*acc(p,Mass,G,N);
                kr1 = dt*v;
                kv2 = dt*acc(p+kr1/2,Mass,G,N);
                kr2 = dt*(v+kv1/2);
                kv3 = dt*acc(p+kr2/2,Mass,G,N);
                kr3 = dt*(v+kv2/2);
                kv4 = dt*acc(p+kr3,Mass,G,N);
                kr4 = dt*(v+kv3);
                p = p + 1/6*(kr1+2*kr2+2*kr3+kr4);
                v = v + 1/6*(kv1+2*kv2+2*kv3+kv4);
            end
            
            %fetch the kinetic and potential energy.
            [kin,pot] = EnergyTracer(p,N,v,Mass,G);
            
            %make a Total kinetic energy vector for plotting.
            E_tot(index) = (kin + pot - E_0) / E_0;
            
            %Calculate the Angular momentum
            L = AngularMomentum(p,N,Mass,v);
            
            %make a angular momentum vector for plotting
            L_t(index) = (L(3)-L_0(3))/L_0(3);
            
            if type~= 2
                %make a momentum vector for plotting (only the norm)
                momentum(:,:,index) = Mass.*v; % momentum of all particles (3xNxtime)
                momentum_norm = vecnorm(nansum(momentum,2),2,1); %(1x1xtime)
                rel_momentum = momentum_norm./vecnorm(momentum(:,end-3,1),2,1); %momentum relative to jupiter
                rel_momentum = permute(rel_momentum,[3,2,1]);
            end
            
            [ecc, semi_m_axis] = eccentricity_sma(p,v,Mass, p);
            ecc = ecc';
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
                    %particle system
                    subplot(2,2,3)
                    plot(p(1,2:end),p(2,2:end),'.k','MarkerSize',20); hold on
                    plot(p(1,1),p(2,1),'*y', 'MarkerSize',20); hold on
                    axis([-1 1 -1 1]*defaultRange*1.1);
                    title(strcat('N =', " ", num2str(sum(Mass~=0)-1)));
                    
                    %plot kuiperbelt if type == 2
                    if type == 2
                        plot(p_k(1,:),p_k(2,:),'.r','MarkerSize',20); hold on
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
        
        xCells{power} = xPower;
        
        saved_var(power) = xPower(1,end);%p(1,2);%vecnorm(v(:,2));%kin + pot;% L(3)-L_0(3);%;
        
        if power>=3
            p2(int_met,power-2) = (saved_var(power-1) - saved_var(power-2))./(saved_var(power)-saved_var(power-1));
        end
    end
    figure(2)
    plot(saved_var)
    hold on
end
legend('1','2','4','5','6','7','8')
%p2(int_methods,:)
disp(log(abs(p2(int_methods,:)))/log(2));
figure(3)
plot(log(abs(p2(int_methods,:)'))/log(2));
legend('1','2','4','5','6','7','8')
hold on
for i = 1:size(int_methods,2)
    text(1,log(abs(p2(int_methods(i),1)'))/log(2),num2str(int_methods(i)))
end
figure(4)
for i = 1:numel(xCells)
    plot(xCells{i,1}(2,:),xCells{i,1}(1,:))
    hold on
end
