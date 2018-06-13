function [chain] = chain_finder(chain, indices)
    for i = chain
        up = find(indices(1,:),i);
        down = find(indices(2,:),i);
        if ismember(chain,up) == 0
            chain = [chain up];
            chain_finder(chain, indices)
        end
        if ismember(chain,down) == 0
            chain = [chain down];
            chain_finder(chain, indices)
        end
    end
end

