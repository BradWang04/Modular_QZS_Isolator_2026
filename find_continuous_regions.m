% 找到连续的QZS区间
function [start_idx, end_idx] = find_continuous_regions(indices)
    if isempty(indices)
        start_idx = [];
        end_idx = [];
        return;
    end
    
    diff_indices = diff(indices);
    breaks = find(diff_indices > 1);
    
    start_idx = [indices(1); indices(breaks + 1)];
    end_idx = [indices(breaks); indices(end)];
end