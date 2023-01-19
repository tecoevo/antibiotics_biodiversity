
function mmod = modifiedmod(a,b)
% Returns rem = remainder of a after dividing by b if rem ~= 0. If rem = 0, returns b. 
% Used to iterate over indices of matrix elements. 
    if mod(a,b) ~= 0
        mmod = mod(a,b);
    elseif mod(a,b) == 0
        mmod = b;
    end
end
