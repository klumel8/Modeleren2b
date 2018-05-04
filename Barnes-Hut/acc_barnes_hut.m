function a = acc_barnes_hut(p, Mass, G, N,theta)
% not ready yet, but currently no time
    %create root of the tree
    basic_tree = tree('0');
    for i = 1:4
        basic_tree = basic_tree.addnode(1,['0',num2str(i)]);
    end
    %define the maximal range
    start_range = 1.1*max(max(abs(p),1),2);
    
    %create trees
    basic_tree = make_tree(basic_tree,2,p,start_range);
    [mass_tree, pos_tree ]= make_mass_pos_tree(basic_tree,p,Mass,start_range);
    
    a = zeros(3,N);
    %per particle, loop through tree, starting from root, and calculate
    %forces
    for i = 1:N
        %maybe works without loop when removing (i)
        distance_tree = pos_tree - p(i);
        distance_norm_tree = distance_tree.treefun(@(x) norm(x,1));
        %calculate which center of masses are far enough
        far_enough = (0.5.^(distance_norm_tree.depthtree-1)).*start_range./distance_norm_tree > theta;
        %calculate the force of each center of mass on this particle
        force_tree = G.*mass_tree.*Mass.*distance_tree./(distance_norm_tree.^3);
        
        % set elements after the first non-zero element to 0?
        
        % sum over first elements from the root which arent 0
        % currently no idea how without loop
        iter = force_tree.depthfirstiterator;
        total_force = force_tree.get(1).*0; % get right size
        for index = iter
            if far_enough(index) && ~far_enough(far_enough.getparent(index))
                %only adds force when parent is not far enough and current
                %node is
                total_force = total_force + force_tree(index);
            end
        end
        
        a(:,i) = total_force./Mass(i);       
    end
    
    
end