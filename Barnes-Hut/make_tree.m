function [new_tree] = make_tree(curr_tree,first_leaf,p,start_range)
%MAKE_TREE create a tree(recursively) based on the positions of particles, starting with
%a square around the origin with size 2*start_range
%input arguments:
%   curr_tree   : (tree) the current tree, on which nodes will be added
%   first_leaf  : (scalar) the first leaf of the tree
%   p           : (3xN) position vector
%   start_range : (scalar) the biggest range to be worked with
%output arguments:
%   new_tree    : (tree) a new tree with one particle in each leaf
%required classes/functions:
%   @tree
    % make startrange a power of 2:
    q = ceil(log(start_range)/log(2));
    start_range = 2^q;
    centers = [.5, -.5,-0.5,0.5,.5, -.5,-0.5,0.5;0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5;0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5];
    new_tree = tree(curr_tree);
    iterator = curr_tree.breadthfirstiterator;
    stop_next = true;
    for i = iterator(first_leaf:end)
        if i ~= 1
            curr_node = curr_tree.get(i);
            indices = str2double(split(curr_node(2:end),''));
            indices = indices(2:end-1);%remove NaN
            n = numel(indices)-1;

            curr_center = sum(0.5.^(0:n).*centers(:,indices),2);
        else
            curr_center = [0;0;0];
            n = 0;
        end
        for j = 1:8
            
            new_center = curr_center+0.5^n.*centers(:,j);
            x_check = [(new_center(1)+centers(1,j)*0.5^n)*start_range;(new_center(1)-centers(1,j)*0.5^n)*start_range];
            x_right = (p(1,:) <= max(x_check) & p(1,:) > min(x_check));
            y_check = [(new_center(2)+centers(2,j)*0.5^n)*start_range;(new_center(2)-centers(2,j)*0.5^n)*start_range];
            y_right = p(2,:) <= max(y_check) & p(2,:) > min(y_check);
            z_check = [(new_center(3)+centers(3,j)*0.5^n)*start_range;(new_center(3)-centers(3,j)*0.5^n)*start_range];
            z_right = p(3,:) <= max(z_check) & p(3,:) > min(z_check);
            all_right = x_right & y_right & z_right;
            
            if any(all_right)
                new_tree = new_tree.addnode(i,[curr_tree.get(i), num2str(j)]);
                stop_next = stop_next & sum(all_right)<=1;
            end
        end
    end
    first_leaf = iterator(end)+1;

    if ~stop_next
        new_tree = make_tree(new_tree,first_leaf,p,start_range);
    end
end