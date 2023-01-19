function product = GenerateAllSequences(n,k) 
%Generates all sequences with entries in alphabet {1,2,...,n} having length k. We do
%this by taking successive cartesian products of {1,2...n} with itself.
%Each sequence is a column in the output `product`. 

currentset = 1:n;

for iteration=1:k-1
    [rowsize,colsize] = size(currentset);
    dummy = zeros(rowsize+1,colsize*n); %dummy will store this iteration's cartesian product
    for i = 1:colsize
        copies_of_col_i = repmat(currentset(:,i), 1, n); %take the first column and generate n copies of it
        products_involving_coli = [copies_of_col_i ; 1:n]; % to each copy, attach a different number from {1,2,...,n}
        for j = 1:n
            dummy(:,j + (i-1)*n) = products_involving_coli(:,j); %add these new copies to dummy
        end    
    end
    currentset = dummy; %take the current value of the process for the next iteration
end

product=currentset;

end
