function [new_tree] = make_tree(curr_tree,first_leaf,p,Mass,start_range)
%MAKE_TREE create a tree(recursively) based on the positions of particles, starting with
%a square around the origin with size 2*start_range
%input arguments:
%   curr_tree   : (tree) the current tree, on which nodes will be added
%   first_leaf  : (scalar) the first leaf of the tree
%   p           : (3xN) position vector
%   Mass        : (1xN) mass vector
%   start_range : (scalar) the biggest range to be worked with
%output arguments:
%   new_tree    : (tree) a new tree with one particle in each leaf
%required classes/functions:
%   @tree
% 
%layout used octants:
%z>0:
% -----------------------
% |    2    |     1     |
% |         |           |
% -----------------------
% |    3    |     4     |
% |         |           |
% -----------------------
%z<= 0:
% -----------------------
% |    6    |     5     |
% |         |           |
% -----------------------
% |    7    |     8     |
% |         |           |
% -----------------------
%example name of node: 015
%0 does not mean anything, used for initiation
%1: first octant (see above)
%5: first octant is divided into another 8 octants, the node represents the
%5th octant(see above) within the first octant( thus 0,0,0 is replaced by
%the center of the first octant)

    %initialize used variables:

    % make startrange a power of 2:
    q = ceil(log(start_range)/log(2));
    start_range = 2^q;
    %define the centers (relative, startrange = 1)
    centers = [ 0.5,-0.5,-0.5, 0.5, 0.5,-0.5,-0.5, 0.5; ...
                0.5, 0.5,-0.5,-0.5, 0.5, 0.5,-0.5,-0.5; ...
                0.5, 0.5, 0.5, 0.5,-0.5,-0.5,-0.5,-0.5];
            
    %initiate tree (copy of curr_tree)
    new_tree = tree(curr_tree);
    %initiate iterator
    iterator = curr_tree.breadthfirstiterator;
    %if stop_next is true at the end of the loop, the recursion will stop
    stop_next = true;
    %loop over leaves, and adds node to a leaf if there is a particle in
    %that leaf
    for i = iterator(first_leaf:end)
        if i ~= 1
            %get the node name, to determine which cell(in positions) it is
            %example name: 0148: first octant, then within that octant the
            %4th octant, then within thát octant the 8 octant
            curr_node = curr_tree.get(i);
            indices = str2double(split(curr_node(2:end),''));
            indices = indices(2:end-1);%remove NaN
            n = numel(indices);

            curr_center = sum(0.5.^(0:(n-1)).*centers(:,indices),2);
            %calculate number of particles in this cell of which the
            %mass is not 0(particles with 0 mass have collided and should
            %left out of the calculation)
            x_check = [(curr_center(1)+0.5^(n));(curr_center(1)-0.5^(n))]*start_range;
            x_right = p(1,:) <= max(x_check) & p(1,:) > min(x_check)& Mass(:)~=0;
            y_check = [(curr_center(2)+0.5^(n));(curr_center(2)-0.5^(n))]*start_range;
            y_right = p(2,:) <= max(y_check) & p(2,:) > min(y_check)& Mass(:)~=0;
            z_check = [(curr_center(3)+0.5^(n));(curr_center(3)-0.5^(n))]*start_range;
            z_right = p(3,:) <= max(z_check) & p(3,:) > min(z_check)& Mass(:)~=0;
            all_right = x_right & y_right & z_right;
            n_particles = sum(all_right);
            
        else
            %for the case the only node is the root
            curr_center = [0;0;0];
            n = 0;
            n_particles = size(p,2) - size(Mass(Mass~=0));
        end
        %only adds nodes if there is more than 1 particle in this cell with
        %mass unequal 0
        if n_particles > 1
            for j = 1:8
                %checks if there is one or more particle with mass ~= 0 in
                %the new node
                new_center = curr_center+0.5^n.*centers(:,j);
                x_check = [(new_center(1)+centers(1,j)*0.5^n);(new_center(1)-centers(1,j)*0.5^n)]*start_range;
                x_right = (p(1,:) <= max(x_check) & p(1,:) > min(x_check))& Mass(:)~=0;
                y_check = [(new_center(2)+centers(2,j)*0.5^n);(new_center(2)-centers(2,j)*0.5^n)]*start_range;
                y_right = p(2,:) <= max(y_check) & p(2,:) > min(y_check)& Mass(:)~=0;
                z_check = [(new_center(3)+centers(3,j)*0.5^n);(new_center(3)-centers(3,j)*0.5^n)]*start_range;
                z_right = p(3,:) <= max(z_check) & p(3,:) > min(z_check)& Mass(:)~=0;
                all_right = x_right & y_right & z_right;
                %only adds node if there is a particle in the next part
                %with mass>0
                if any(all_right)  
                    %add node:
                    new_tree = new_tree.addnode(i,[curr_tree.get(i), num2str(j)]);
                    
                    %stops the iteration if all new nodes have only 1 particle
                    stop_next = stop_next & sum(all_right)==1; 
                end
            end
        end
    end
    %gets next first leaf
    first_leaf = iterator(end)+1;
    
    %stop when each particle has a cell 
    if ~stop_next 
        new_tree = make_tree(new_tree,first_leaf,p,Mass,start_range);
    end
    disp(new_tree.tostring)

end