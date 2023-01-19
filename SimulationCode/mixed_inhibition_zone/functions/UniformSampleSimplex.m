function array = UniformSampleSimplex(n,k)
% Generates k points uniformly sampled from the (n-1)simplex in R^n i.e., 
% points whose coordinates all sum to 1. Returns n x k matrix, 
% each column is a different point. 
% method taken from Donald B. Rubin. "The Bayesian Bootstrap." Ann. Statist. 9 (1) 130 - 134, January, 1981. https://doi.org/10.1214/aos/1176345338

array = zeros(n,k);
for i = 1:k
    initialize = sort([rand(1,n-1), 0, 1]);
    draw = zeros(1,n);    
    for coord = 1:n
        draw(coord) = initialize(coord+1) - initialize(coord);
    end
    array(:,i) = draw;
end

end