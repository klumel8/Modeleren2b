function a = acc_barnes_hut(p, Mass, G,N,theta)
% does not always work, maybe because of collisions

    %create root of the tree
    basic_tree = tree(0);

    %define the maximal range
    start_range = 1.1*max(max(max(abs(p))));
    %start_range = 8e11;
    % make startrange a power of 2:(maybe is not needed)
    q = ceil(log(start_range)/log(2));
    start_range = 2^q;

    %create trees
    basic_tree = make_tree(basic_tree,1,p,Mass,start_range);
    [mass_tree, pos_tree ]= make_mass_pos_tree(basic_tree,p,Mass,start_range);
    
    %calculate distance between particles and center of masses
    distance_tree = pos_tree - p;
    distance_norm_tree = distance_tree.treefun(@(x) vecnorm(x));
    
    %calculate which center of masses are far enough, s/d < theta: s =
    %longest distance in one cell, d = distance between center of mass of
    %the cell and the particle
    depthNodes = cellfun('length',basic_tree.Node)-1;
    tempNode = 0.5.^(depthNodes-1);
    tempTree = basic_tree;
    tempTree.Node = num2cell(tempNode);
    
    tempNode = ((3.*((0.5.^(depthNodes-1)).*start_range).^2).^0.5)./cell2mat(distance_norm_tree.Node);
    far_enough_Node = tempNode < theta;
    %oldTemp = 0.5.^(distance_norm_tree.depthtree-1);
    far_enough_Node = mat2cell(far_enough_Node,ones(1,size(far_enough_Node,1)),size(far_enough_Node,2));
    far_enough = basic_tree;
    far_enough.Node = far_enough_Node;
    
    
    %far_enough = (3.*((tempTree).*start_range).^2).^(0.5)./distance_norm_tree < theta;
%     far_enough = ((0.5.^(distance_norm_tree.depthtree-1)).*start_range)./distance_norm_tree < theta;
    
    
    depth_iter = 1:numel(far_enough_Node);
    %set leaves to 1, to ensure that there is a force to calculate
    %set far_enough to 0 if one or more of its children are 0
    for i = depth_iter(end:-1:1)
        if far_enough.isleaf(i)
            %get the index of the particle in this leaf
            [~, min_index] = min(distance_norm_tree.get(i));
            
            far_enough = far_enough.set(i,(1:N ~= min_index)); % 0 at the min index, 1 elsewhere
        else
            for child = far_enough.getchildren(i)
                far_enough = far_enough.set(i,far_enough.Node{i}.*far_enough.Node{child});
            end
%         far_enough = far_enough.set(far_enough.getparent(i),far_enough.get(far_enough.getparent(i)).*far_enough.get(i));
            if any(far_enough.Node{i})
                %disp(far_enough.Node{i})
            end
        end
    end

    
    %set root to 0, to ensure there some nodes are used to calculate the
    %force
    far_enough = far_enough.set(1,zeros(size(Mass)));
    
    %calculate the force of each center of mass on this particle
    %sasha: calculate force of CoM's on each other if far enough away, then
    %distribute force over particles
    %nando: per particle: calculate force of CoM's which are far enough
    %away
    %nando: dont know if it is possible without forloops
    
    
    %goes wrong: (i think its only this part, i could be wrong)
    a = zeros(3,N);
    a_T = zeros(3,N);
    iter = 1:numel(far_enough.Node);
    
    %loop over CoM
    for k = iter(2:end) %exclude root from iteration, because it doesnt have a parent and should never be used to calculate the force
        far_enough.getparent(k);
        
        parent_node = far_enough.get(far_enough.getparent(k));
        node = far_enough.get(k); %node = 1xN matrix with 1 at place i when this CoM is far enough from particle i
        
        %{
        for i = 1:N %loop over particles
            %if the CoM is far enough, and if parent is not far enough, calculate the acceleration
            
            if node(i) & ~parent_node(i)
                %direction of acceleration is right
                %maybe calculate distance before forloop?
                distance = pos_tree.get(k) - p(:,i); %vector in the direction of the CoM
                a(:,i) = a(:,i) + G.*mass_tree.get(k).*distance./(vecnorm(distance).^3);
            end
                
        end
        %}
        distance = (pos_tree.get(k)-p(:,:)).*repmat(node.*(~parent_node),3,1);
        acc_added = G.*mass_tree.get(k).*distance./(vecnorm(distance,2,1).^3);
        acc_added(isnan(acc_added)) = 0;
        a = a + acc_added;
    end
    
        
    %sasha: 
    %suppose we calculate the force F from CoM a1 to CoM b1, and a2 is a
    %particle (or CoM) contributing to a1:
    % then we have vectors b1-a1, and a2-a1 -> b1 -a1 - (a2- a1) = b1-a2 =
    % direction of the force from b1 acting on a2. the size of the force =
    % F*m_a1/m_a2 -> F_a2 = |F|*m_a1/m_a2 * (b1-a2)/|b1-a2|

    

    
    
end