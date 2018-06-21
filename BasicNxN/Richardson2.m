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
int_methods = [7];

for int_met = int_methods
    start_dt = 3600*24*7*52*10;
    T = start_dt;
    saved_var = [];
    P = 5;
    xCells = cell(P,1);
    for power = 1:P
        power
        dt = start_dt*2^(-(power-2));
        %Particles in our model;
        % 1 = early solar system
        % 2 = solar system and Kuyper belt
        % 3 = sphere
        type = 2;

        make_movie = false;
        cycle_count = 0;
        d_theta_old = 0;
        
        
        barnes_hut = false;
        theta = 0.5;%0 to test acc calculation: all particles are indiviually used,
        %should be the same as without barnes hut
        
        % universal parameters
        G = 6.67408*10^-11; % [Nm^2kg^-2]
        AU = 1.49597871e11; % [m]
        
        rng(121) %rng(seed): Used to control random number generation
        
        if type == 2 % solar system and Kuyper belt
            defaultRange = 50*AU; % [m]
            N = 1e4; % Dummy variable
            N_k = 0; % particles in kuiper belt
            %dt = 3600*24*7*52*3; % in seconds
            T = start_dt*2^12; % in seconds
            [Mass, p, v, N] = initialConditions(defaultRange,N,2);
            [p_k, v_k, Mass_k] = kuiperbelt(N_k, p);
            max_orbit_length = 20000; %determines how much of the orbit of a single particle is shown
            particle = 2; % plot path of Pluto
            
            kuipercollisions = false;
        end
        
        
        % plotting configuration
        plot_system = false;     %plot the particle system
        plot_ecc_a = false;      %plot eccentricity vs semi major axis
        plot_ang_mom = false;    %plot the angular momentum
        plot_momentum = false;   %plot the momentum, relative to jupiter(only for type ==2)
        plot_RV = false;          %plot the range vs the speed
        plotting = false;        %plot anything at all
        plot_hist = false;
        
        fps = 1/3;
        TstepsPframe = 1/4;
        frames = floor(T/(TstepsPframe*dt))+1;
        if make_movie
            F(frames) = struct('cdata',[],'colormap',[]);
        end
                
        
        %remove the particles every [remove_index] timesteps
        remove_index = 10;
        removing = true;
        
        %define acc function:
        if barnes_hut
            acc_fun = @(p,Mass,N) acc_barnes_hut(p,Mass,G,N,theta);
        else
            acc_fun = @(p,Mass,N) acc(p,Mass,G,N);
        end
        
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
        
        single_p = [];
        % a timer so we dont plot too often and slow down the script
        tic;
        for t = 0:dt:T
            xPower(:,t/dt+1) = [(p(1,2)^2+p(2,2)^2)^0.5;t];
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
            %remove_crit = vecnorm(p)>100*defaultRange;
            remove_crit = sqrt(sum(p.^2,1))>100*defaultRange;
            
            if any(remove_crit)
                disp('too far')
                indices = find(remove_crit)
                disp(['Number of removed particles: ',num2str(numel(indices))])
                p(:,indices) = p(:,indices)./20;%bring particle back in the system
                v(:,indices) = repmat([0;0;0],1,numel(indices));%set velocity to 0
                Mass(indices) = 0;%set mass to 0
            end
            
            %read fo.m first, but keeps track of whether there was a collision.
            c = col(p,v,Mass,N,dt);
            %#BUG will crash if multiple collisions in one timestep
            %check if the collision vector is empty
            if max(max(c)) > 0
                disp('collision')
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
            
            if type == 2 && kuipercollisions
                c_k = col_kuiper(p,p_k,v,v_k,Mass,dt);
                if any(any(c_k))
                    disp('collision(kuiperbelt particle)')
                    %find indices of collided particles
                    %re-rank the collision indexes
                    indices = [mod(find(c_k),N_k)'; ceil(find(c_k)/N)'];
                    indices(1,:) = (indices(1,:)==0)*N_k + indices(1,:);
                    removing_indices = zeros(1,N_k);
                    removing_indices(indices(1,:)) = 1;
                    %remove kuiper belt particles which collide
                    p_k = p_k(:,~removing_indices);
                    v_k = v_k(:,~removing_indices);
                    if t>0
                        clearvars plot_p_k
                    end
                    if int_met == 7
                        a_k = a_k(:,~removing_indices);
                    end
                    N_k = size(p_k,2);
                end
            end
            
            
            if int_met == 1
                %first order (newton forward)
                a = acc_fun(p,Mass,N);
                p = p + v*dt + a * dt^2 / 2;
                v = v + a * dt;
            elseif int_met == 2
                %Runge Kutta 2 conserves angular momentum?
                k1 = dt^2*acc_fun(p + (1/2)*dt*v,Mass,N);
                p = p + v*dt + k1;
                v = v + k1/dt;
            elseif int_met == 4
                %Runge Kutta 4
                if type ~= 2
                    k1 = dt^2*acc_fun(p,Mass,N);
                    k2 = dt^2*acc_fun(p + 0.5*dt*v + 1/8*k1,Mass,N);
                    k3 = dt^2*acc_fun(p + dt*v + .5*k2,Mass,N);
                    
                    p = p + v*dt + 1/6*(k1+2*k2);
                    v = v + 1/(6*dt)*(k1+4*k2+k3);
                end
                %a_k = kuiperacc(p,p_k,Mass);
                if type == 2
                    [k1_k,k1_pl] = kuiperacc(p,p_k,Mass,Mass_k);
                    k1_k = dt^2*k1_k;
                    k1 = dt^2*acc_fun(p,Mass,N);
                    k1 = k1+dt^2*k1_pl;
                    
                    [k2_k,k2_pl] = kuiperacc(p + 0.5*dt*v + 1/8*k1,p_k + 0.5*dt*v_k + 1/8*k1_k ,Mass,Mass_k);
                    k2_k = dt^2*k2_k;
                    k2 = dt^2*acc_fun(p + 0.5*dt*v + 1/8*k1,Mass,N);
                    k2 = k2+dt^2*k2_pl;
                    
                    [k3_k,k3_pl] = kuiperacc(p + dt*v + .5*k2,p_k + dt*v_k + .5*k2_k,Mass,Mass_k);
                    k3_k = dt^2*k3_k;
                    k3 = dt^2*acc_fun(p + dt*v + .5*k2,Mass,N);
                    k3 = k3+dt^2*k3_pl;
                    
                    p_k = p_k + v_k*dt + 1/6*(k1_k+2*k2_k);
                    v_k = v_k + 1/(6*dt)*(k1_k+4*k2_k+k3_k);
                    
                    p = p + v*dt + 1/6*(k1+2*k2);
                    v = v + 1/(6*dt)*(k1+4*k2+k3);
                end
                
                
            elseif int_met == 5
                %Runge Kutta 5a see file I (floris) send over whatsapp
                k1 = dt^2*acc_fun(p,Mass,N);
                k2 = dt^2*acc_fun(p + (1/4)*dt*v + (1/32)*k1,Mass,N);
                k3 = dt^2*acc_fun(p + (7/10)*dt*v - (7/1000)*k1 + (63/250)*k2,Mass,N);
                k4 = dt^2*acc_fun(p + dt*v + (2/7)*k1 + (3/14)*k3,Mass,N);
                p = p + dt*v + (1/14)*k1 + (8/27)*k2 + (25/189)*k4;
                v = v + (1/dt)*((1/14)*k1 + (32/81)*k2 + (250/567)*k3 + (5/54)*k4);
            elseif int_met == 6
                %Runge Kutta 6.
                k1 = dt^2*acc_fun(p,Mass,N);
                k2 = dt^2*acc_fun(p + 1/4*dt*v + 1/32*k1,Mass,N);
                k3 = dt^2*acc_fun(p + 1/2*dt*v - k1/24 + k2/6,Mass,N);
                k4 = dt^2*acc_fun(p + 3/4*dt*v + k1*3/32 + k2/8 + k3/16,Mass,N);
                k5 = dt^2*acc_fun(p + 3/7*dt*v - k1/14 + k3/7,Mass,N);
                p = p + dt*v + (7*k1 +24*k2 + 6*k3 + 8*k4)/90;
                v = v + (7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5)/(90*dt);
            elseif int_met == 7
                if type ~=2
                    if t == 0
                        %initialize acceleration for leapfrog
                        a = acc_fun(p,Mass,N);
                    end
                    %leapfrog
                    v = v + a*dt/2;
                    p = p + dt*v;
                    a = acc_fun(p,Mass,N);
                    
                    v = v + a*dt/2;
                end
                
                if type == 2
                    if t == 0
                        %initialize acceleration for leapfrog
                        [a_k,a_pl] = kuiperacc(p,p_k,Mass,Mass_k);
                        a = acc_fun(p,Mass,N)+a_pl;
                    end
                    %leapfrog
                    v = v + a*dt/2;
                    p = p + dt*v;
                    a = acc_fun(p,Mass,N);
                    v = v + a*dt/2;
                    
                    v_k = v_k + dt/2*a_k;
                    p_k = p_k + dt*v_k;
                    [a_k,a_pl] = kuiperacc(p,p_k,Mass,Mass_k);
                    a = a + a_pl;
                    v_k = v_k + a_k*dt/2;
                    
                end
            end
            
            
            %fetch the kinetic and potential energy.
            [kin,pot] = EnergyTracer(p,N,v,Mass,G);
            
            %make a Total kinetic energy vector for plotting.
            %     E_tot(index) = (kin + pot - E_0) / E_0;
            
            %Calculate the Angular momentum
            L = AngularMomentum(p,N,Mass,v);
            
            %make a angular momentum vector for plotting
            %     L_t(index) = (L(3)-L_0(3))/L_0(3);
            
            %     if type == 2
            %         %make a momentum vector for plotting (only the norm)
            %         momentum(:,:,index) = Mass.*v; % momentum of all particles (3xNxtime)
            %         momentum_norm = vecnorm(nansum(momentum,2),2,1); %(1x1xtime)
            %         rel_momentum = momentum_norm./vecnorm(momentum(:,end-3,1),2,1); %momentum relative to jupiter
            % %         rel_momentum = momentum_norm./vecnorm(momentum(:,6,1),2,1); %momentum relative to jupiter
            %         rel_momentum = permute(rel_momentum,[3,2,1]);
            %     end
            
            [ecc, semi_m_axis] = eccentricity_sma(p,v,Mass, p);
            
            
            ecc = ecc';
            
            semi_m_axis = semi_m_axis';
            if type == 2
                [ecc_kuiper, semi_m_axis_kuiper] = eccentricity_sma(p_k,v_k,Mass,p);
                ecc_kuiper = ecc_kuiper';
                
                semi_m_axis_kuiper = semi_m_axis_kuiper';
            end
            d_theta = atan(p(2,2)/p(1,2));
            d_theta = d_theta - pi*(p(1,2)<0)+pi/2;
            
            A = [cos(d_theta), sin(d_theta); -sin(d_theta), cos(d_theta) ];
            single_p = [single_p, A*p(1:2,particle)];
            
            if size(single_p,2)>max_orbit_length
                single_p = single_p(:,2:end);
            end
            
            
            
            
            
            
            %when plotting too often this can drastically slow down the script. Plotting once every 200 timesteps help speeding this up IFF the plotting is bottlenecking the script
            %only plot when plotting = true, (saves time)

            curr_time = toc;
            if (plotting && (mod(t,TstepsPframe*dt)==0)  && (curr_time>1/fps || t == 0))
                
                if t == 0
                    figure(1)
                end
                
                if ~plot_ang_mom
                    figure(2);
                    %Produce figures with a LaTeX interpreter
                    set(0,'defaulttextinterpreter','latex');
                    set(0,'defaultaxesfontsize',14);
                    set(gca, 'ticklabelinterpreter','latex');
                    plot_p(1:2,:) = A*p(1:2,:);
                    
                    ax_single = gca;
                    plot(single_p(1,:), single_p(2,:),'-b','LineWidth',0.05);
                    ax_single.NextPlot = 'add'; %Hold on, maar dan dat de assen ook bewaren
                    
                    
                    plot(plot_p(1,2:end),plot_p(2,2:end),'.k','MarkerSize',20); hold on;
                    plot(plot_p(1,1),plot_p(2,1),'*y', 'MarkerSize',20);
                    axis([-1.2 1.2 -1.2 1.2]*defaultRange);
                    ax_single.NextPlot = 'replaceChildren'; %Houdt dezelfde assen nu ook bij vervolgplots
                    title({'Resonances in Kuiper belt',strcat('time: ',num2str(round(t/31556926,1)),' y')})
                    xlabel('$x$ [m]'); ylabel('$y$ [m]')
                    
                end
                
                drawnow
                
                if make_movie
                    F(t/(TstepsPframe*dt)+1) = getframe(gcf);curr_tree.Node{i};
                end
                tic;
            end
        end
        xCells{power} = xPower;
        
        saved_var(power) = xPower(1,end);%p(1,2);%vecnorm(v(:,2));%kin + pot;% L(3)-L_0(3);%;
        
        if power>=3
            p2(int_met,power-2) = (saved_var(power-1) - saved_var(power-2))./(saved_var(power)-saved_var(power-1));
        end
    end
    figure(3)
    plot(saved_var)
    hold on
end

order = log(abs(p2(int_methods,:)))/log(2);
disp(order);
figure(5)
plot(order');
legend('4')
hold on
for i = 1:size(int_methods,2)
    text(1,log(abs(p2(int_methods(i),1)'))/log(2),num2str(int_methods(i)))
end
figure(4)
for i = 1:numel(xCells)
    plot(xCells{i,1}(2,:),xCells{i,1}(1,:))
    hold on
end
