function array = struct2array(struct)
% converts structure arrays containing sym numbers to an array containing
% the corresponding numbers in double precision. 

number_of_fixedpts = numel(struct.X1);  % vpasolve may return more than one fixed point
array = zeros(number_of_fixedpts,1);         %initialize the output array

fn = fieldnames(struct);    %get all fieldnames of the structure array
for k=1:numel(fn)  
    array(:,k) = double(struct.(fn{k}));
end
array = array';

end