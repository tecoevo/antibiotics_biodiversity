function finalprobability = InclusionExclusion(iterable)
% iterable must be a cell array
% Given a collection of probabilities of possibly overlapping events 
% (inputted as a 1-dimensional iterable object) this function computes the probability 
% of the union of these events using the principle of inclusion-exclusion. 

finalprobability = 0;

len_iter = numel(iterable); % number of events
indicator_array = 1:len_iter; 

for subset_size = indicator_array
    contribution_of_size = 0; % contribution of subsets of this size i.e., subset_size
    subsets = nchoosek(indicator_array,subset_size)'; %columns are different subsets
    for subset = subsets
        contribution_of_subset = 1;
        for index = subset' 
            contribution_of_subset = contribution_of_subset*iterable{index};
        end
        contribution_of_size = contribution_of_size + contribution_of_subset;
    end
    finalprobability = finalprobability + (-1)^(subset_size - 1)*contribution_of_size;
end

end