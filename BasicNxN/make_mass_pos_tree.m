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
    
    %define the centers(relative, startrange = 1)
    centers = [ 0.5,-0.5,-0.5, 0.5, 0.5,-0.5,-0.5, 0.5; ...
                0.5, 0.5,-0.5,-0.5, 0.5, 0.5,-0.5,-0.5; ...
                0.5, 0.5, 0.5, 0.5,-0.5,-0.5,-0.5,-0.5];
    %initiate the trees(copy of basic_tree, with either 0 as data, or [0;0;0]
    mass_tree = tree(basic_tree,0);
    pos_tree = tree(basic_tree,[0;0;0]);
    %initiate iterator
    iterator = 1:numel(basic_tree.Node);
    
    %fill leaves in mass_tree and pos_tree
    for i = iterator(end:-1:1)
        if basic_tree.isleaf(i)
            %get node name
            curr_node = basic_tree.Node{i};
            %temp = split(curr_node(2:end),'');
            %indices = str2double(temp);
            %indices = indices(2:end-1);%remove NaN
            indices = curr_node(2:end);
            n = numel(indices);
            %get center position
            curr_center = sum(0.5.^(0:(n-1)).*centers(:,indices),2);

            %get the particle which is in this node(neglecting the
            %particles with Mass == 0
            x_check = [(curr_center(1)+0.5^(n));(curr_center(1)-0.5^(n))]*start_range;
            x_right = p(1,:) <= max(x_check) & p(1,:) > min(x_check);
            y_check = [(curr_center(2)+0.5^(n));(curr_center(2)-0.5^(n))]*start_range;
            y_right = p(2,:) <= max(y_check) & p(2,:) > min(y_check);
            z_check = [(curr_center(3)+0.5^(n));(curr_center(3)-0.5^(n))]*start_range;
            z_right = p(3,:) <= max(z_check) & p(3,:) > min(z_check);
            all_right = x_right & y_right & z_right  & Mass~=0;
            
            %sum is only needed if we are going to put multiple particles
            %in one cell: (maybe that wont work properly, it isnt
            %implemented yet)
            mass_tree.Node{i} = sum(Mass(all_right)); 
            %sums are only needed if we are going to put multiple particles
            %in one cell, otherwise only p(:,all_right) is needed:(see
            %comment above)
            if any(p(3,all_right)~=0)
                %disp(p(:,all_right))
            end
            pos_tree.Node{i} = sum(Mass(all_right).*p(:,all_right),2)./sum(Mass(all_right));
        end 
    end
    %fill mass_tree by summing over the nodes
    mass_tree = mass_tree.recursivecumfun(@(x) sum(x,2));
    %make a tree with m*r, to easily get the center of mass

    mass_pos_tree = pos_tree.treefun2(mass_tree,@(r,m) m.*r);
    mass_pos_tree = mass_pos_tree.recursivecumfun(@(x) sum(x,2),3);

    %fill the pos_tree (center of masses)
    pos_tree = mass_tree.treefun2(mass_pos_tree,@(m_tot,mr) mr./m_tot);
end
