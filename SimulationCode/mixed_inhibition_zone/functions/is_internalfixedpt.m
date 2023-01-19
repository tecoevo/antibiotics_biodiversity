function logical = is_internalfixedpt(point)
% takes fixed points only as column vectors
% first check if all coordinates are non-zero - we are looking for
% all-species equilibrsym(ia
point = point';

if prod(point)~=0
    interior = 1;
elseif prod(point)==0
    interior = 0;
end

% now check if all coordinates of point are in [0,1]. If not, they cannot
% be abundances.
abundance = 1;
for coordinate = point
    if 0<coordinate && coordinate<1
        abundance = abundance*1;
    elseif coordinate<0 || coordinate>1
        abundance = abundance*0;
    end
end
logical = interior*abundance; % 0 if the point is either not an abundance or not an internal fixed point

end
