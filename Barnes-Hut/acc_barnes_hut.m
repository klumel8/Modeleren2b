function a = acc_barnes_hut(p, Mass, G,N,theta)
% does not always work, maybe because of collisions

    %create root of the tree
    basic_tree = tree('0');

    %define the maximal range
    start_range = 1.1*max(max(max(abs(p))));

    %create trees
    basic_tree = make_tree(basic_tree,1,p,Mass,start_range);
    [mass_tree, pos_tree ]= make_mass_pos_tree(basic_tree,p,Mass,start_range);
    
    %calculate distance between particles and center of masses
    distance_tree = pos_tree - p;
    distance_norm_tree = distance_tree.treefun(@(x) norm(x,1));
    %calculate which center of masses are far enough, s/d < theta: s =
    %longest distance in one cell, d = distance between center of mass of
    %the cell and the particle
    far_enough = (3.*((0.5.^(distance_norm_tree.depthtree-1)).*start_range).^2).^(0.5)./distance_norm_tree < theta;
%     far_enough = ((0.5.^(distance_norm_tree.depthtree-1)).*start_range)./distance_norm_tree < theta;
    %setting leaves to 1, to ensure that there is a force to calculate:
    for i = far_enough.breadthfirstiterator
        if far_enough.isleaf(i)
            far_enough = far_enough.set(i,ones(size(Mass))==1);
        end
    end
    %set root to 0, to ensure there some nodes are used to calculate the
    %force
    far_enough = far_enough.set(1,ones(size(Mass))==0);
    %calculate the force of each center of mass on this particle
    %sasha: calculate force of CoM's on each other if far enough away, then
    %distribute force over particles
    %nando: per particle: calculate force of CoM's which are far enough
    %away
    %nando: dont know if it is possible without forloops
    depth_tree = far_enough.depthtree;
    a = zeros(3,N);
    for i = 1:N
        depth = 0;
        skip_next = false;

        iter = far_enough.depthfirstiterator;
        for k = iter
            if depth >= depth_tree.get(k)
                %set skip_next to false when back at the 'top'(lower depth
                %than previous used node
                skip_next = false;
            end
            if far_enough.isleaf(k)
                if skip_next
                    skip_next = false;
                    continue
                end
            end
            if skip_next
                continue
            end
            node = far_enough.get(k);
            %if the CoM is far enough, calculate the acceleration
            if node(i)
                %skip next node until at the same depth again
                skip_next = true;
                depth = depth_tree.get(k);
                %direction of acceleration is right
                %maybe calculate before forloop?
                distance = pos_tree(k) - p(i);
                a(i) = a(i) + G.*mass_tree.get(k).*distance./(norm(distance).^3);
            end
        end
    end
    
        
    %old:
%     force_tree = G.*mass_tree.*Mass.*distance_tree./(distance_norm_tree.^3);
% 
%        
%     % sum over first elements from the root which arent 0
%     % currently no idea how without loop
%     %dont know if it works properly...
%     
%     iter = force_tree.depthfirstiterator;
%     total_force = force_tree.get(1).*0; % get right size
%     for index = iter(2:end)
%         %only adds force when parent is not far enough and current
%         %node is
%         add_force = (far_enough.get(index) & ~far_enough.get(far_enough.getparent(index)));
%         total_force = total_force + force_tree.get(index).*add_force;
%     end
%         
%     a = total_force./Mass;   

    
    
end