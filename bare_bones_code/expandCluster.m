function [k, point] = expandCluster1(ii, point, k) % pointIndex, clusterID, pointStructure, neighborNumber

if nargin < 3
    k = 3;
end

if k < (point{1}.rLimit)  % Until recursion limit is exceeded
    k = k + 1;
    for j = 1:point{ii}.neighbor.n
        old = point{point{ii}.neighbor.I(j)}.ID;
        point{point{ii}.neighbor.I(j)}.ID = point{ii}.ID;
        if old == 0
            [k, point] = expandCluster1(point{ii}.neighbor.I(j),point, k);
        else
            if (old > 0) && (old ~= point{ii}.ID) % If neighbor is in an older cluster
                gold = point{ii}.ID;
                for p = 1:point{1}.N
                    if point{p}.ID == gold
                        point{p}.ID = old;
                    end
                end
                if point{point{ii}.neighbor.I(j)}.neighbor.n < point{1}.NN
                    H = 'border point';
                else
                    H = 'maximum recursion point';
                end
                G = sprintf('Exception: cluster %d is being merged with cluster %d through %s', gold, old, H);
                disp(G);
            end
        end
    end
else
   point{ii}.ID = 0; 
end

return
