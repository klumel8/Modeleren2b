function cumtree = recursivecumfun(obj, fun,size)
%% RECURSIVECUMFUN  Create a tree where content is calculated recursively from children.
%
% cumtree = RECURSIVECUMFUN(obj, fun) generate a new tree by applying
% recursively at each node the function 'fun' to the content of the
% children nodes. The traversal goes from bottom to top, from leaves to
% root. Keeping this in mind, one understands then that the new tree
% content only depends on the content of the leaves of the source tree.
% 
% EXAMPLE:
%
% extree = tree.example;
% dt = extree.depthtree;
% dt.tostring
% ct = dt.recursivecumfun(@sum);
% ct.tostring
    if nargin == 2
        size = 1;
    end

    % Prepare a blank tree for holding values
    cumtree = tree(obj);
    descend(1);
    
    
    
    
    function val = descend(n) % current node index
        
        if cumtree.isleaf(n)
            val = cumtree.Node{n};
            
        else
           
            children = cumtree.getchildren(n);
            nChildren = numel(children);
            stock = NaN(size,nChildren);
            
            for i = 1 : nChildren
               
                child = children(i);

                stock(:,i) = descend(child);

            end
            val = fun(stock);
            % Store value in new tree
            cumtree.Node(n) = { val };
        end

    end



end