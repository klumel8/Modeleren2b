function [mass_tree, pos_tree] = make_mass_pos_tree(basic_tree,p,Mass,start_range)
%MAKE_MASS_POS_TREE make a tree with total mass of each cell and a tree
%with the center of mass of each cell
%input arguments:
%   basic_tree  : (tree) the used tree structure
%   p           : (3xN) position vector
%   Mass        : (1xN) mass vector
%   start_range  : (scalar) the range used to make the basic_tree
%output arguments:
%   mass_tree   : (tree) a tree with the total mass of each cell
%   pos_tree    : (tree) a tree with the center of mass of each cell
%required classes/functions:
%   @tree
    %make start_range a power of 2
    q = ceil(log(start_range)/log(2));
    start_range = 2^q;
    
    centers = [ 0.5,-0.5,-0.5, 0.5, 0.5,-0.5,-0.5, 0.5; ...
                0.5, 0.5,-0.5,-0.5, 0.5, 0.5,-0.5,-0.5; ...
                0.5, 0.5, 0.5, 0.5,-0.5,-0.5,-0.5,-0.5];
    mass_tree = tree(basic_tree,0);
    pos_tree = tree(basic_tree,[0;0;0]);
    iterator = basic_tree.breadthfirstiterator;
    %doesnt fill all leaves, and sometimes puts 2 particles in 1 leaf
    for i = iterator(end:-1:1)
        if basic_tree.isleaf(i)
            curr_node = basic_tree.get(i);
            indices = str2double(split(curr_node(2:end),''));
            indices = indices(2:end-1);%remove NaN
            n = numel(indices)-1;
            curr_center = sum(0.5.^(0:n).*centers(:,indices),2);

            x_check = [(curr_center(1)+0.5^(n+1));(curr_center(1)-0.5^(n+1))]*start_range;
            x_right = p(1,:) <= max(x_check) & p(1,:) > min(x_check);
            y_check = [(curr_center(2)+0.5^(n+1));(curr_center(2)-0.5^(n+1))]*start_range;
            y_right = p(2,:) <= max(y_check) & p(2,:) > min(y_check);
            z_check = [(curr_center(3)+0.5^(n+1));(curr_center(3)-0.5^(n+1))]*start_range;
            z_right = p(3,:) <= max(z_check) & p(3,:) > min(z_check);
            all_right = x_right & y_right & z_right;

            mass_tree = mass_tree.set(i,Mass(all_right));
            pos_tree = pos_tree.set(i,p(:,all_right));
        end 
    end
    mass_tree = mass_tree.recursivecumfun(@(x) sum(x,2));
    mass_pos_tree = pos_tree.treefun2(mass_tree,@(r,m) m.*r);
    mass_pos_tree = mass_pos_tree.recursivecumfun(@(x) sum(x,2),3);
    pos_tree = mass_tree.treefun2(mass_pos_tree,@(m_tot,mr) mr./m_tot);
end
