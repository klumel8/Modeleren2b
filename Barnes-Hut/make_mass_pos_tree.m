function [mass_tree, pos_tree] = make_mass_pos_tree(basic_tree,p,Mass,startrange)
%MAKE_MASS_POS_TREE make a tree with total mass of each cell and a tree
%with the center of mass of each cell
%input arguments:
%   basic_tree  : (tree) the used tree structure
%   p           : (3xN) position vector
%   Mass        : (1xN) mass vector
%   startrange  : (scalar) the range used to make the basic_tree
%output arguments:
%   mass_tree   : (tree) a tree with the total mass of each cell
%   pos_tree    : (tree) a tree with the center of mass of each cell
%required classes/functions:
%   @tree
    centers = [.5, -.5,-0.5,0.5;0.5,0.5,-0.5,-0.5;0,0,0,0];
    mass_tree = tree(basic_tree,0);
    pos_tree = tree(basic_tree,[0;0;0]);
    iterator = basic_tree.breadthfirstiterator;
    disp((iterator))
    for i = iterator(end:-1:1)
        if basic_tree.isleaf(i)
            curr_node = basic_tree.get(i);
            indices = str2double(split(curr_node(2:end),''));
            indices = indices(2:end-1);%remove NaN
            n = numel(indices);

            curr_center = sum(0.5.^(0:n-1).*centers(:,indices),2);
            x_check = [(curr_center(1)+0.5^n)*startrange,(curr_center(1)-0.5^n)*startrange];
            x_right = (p(1,:) <= max(x_check) & p(1,:) > min(x_check));
            y_check = [(curr_center(2)+0.5^n)*startrange,(curr_center(2)-0.5^n)*startrange];
            y_right = p(2,:) <= max(y_check) & p(2,:) > min(y_check);
            z_check = [(curr_center(3)+0.5^n)*start_range,(curr_center(3)-0.5^n)*start_range];
            z_right = p(3,:) <= max(z_check) & p(3,:) > min(z_check);
            all_right = x_right & y_right & z_right;
            mass_tree = mass_tree.set(i,Mass(all_right));
            pos_tree = pos_tree.set(i,p(:,all_right));
        end 
    end
    mass_tree = mass_tree.recursivecumfun(@(x) sum(x));
    mass_pos_tree = mass_tree.treefun2(pos_tree,@(m,r) m*r);
    mass_pos_tree = mass_pos_tree.recursivecumfun(@(x) sum(x));
    pos_tree = mass_tree.treefun2(mass_pos_tree,@(m_tot,mr) mr/m_tot);
end
