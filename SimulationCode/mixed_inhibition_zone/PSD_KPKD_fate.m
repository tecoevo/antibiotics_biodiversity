N = 3; % number of species
M = 1; % number of antibiotics
g = repelem(1,N); % base fitness of each species
timesteps = 1000;

% metabolic costs of:
c_res = repelem(0.05, M);  %intrinsic resistance
c_prod = 3*c_res; %production
c_deg = 2.1*c_res;  %resistance via degradation

max_parameter_value = 2;
step_size = 0.5;
KP_array = step_size:step_size:max_parameter_value;
KD_array = step_size:step_size:max_parameter_value;
array_length = max_parameter_value/step_size;
result = zeros(array_length);

P = [1;0;0];
S = [0;1;0];
D = [0;0;1];
R = [0;0;0];

fixedpts = zeros(N,1);
i=1; %to increase the size of the array fixedpts

colors = {'red', 'green', 'blue'};
labels = {'P', 'S', 'D'};

samples_of_simplex = 300;

for i = 1:length(KP_array)
    for j = 1:length(KD_array)
        X = sym('X', [1 N], 'real');
        modeleqn = symmodel(X,N,M,P,S,D,R,KD_array(j),KP_array(i),g,c_prod,c_deg,c_res)-X;
        eqns = modeleqn == repelem(0,N);   % equations to solve
        ends = 0;
        for iter = 1:samples_of_simplex
            init = UniformSampleSimplex(N,1);
            trajectory = TimeEvolution(init, N,M,P,S,D,R,KD_array(j),KP_array(i),g,c_prod,c_deg,c_res, timesteps);
            final_pt = trajectory(:,timesteps); 
            if isequal([1;0;0], round(final_pt))
                ends = ends + 1;
            end
            if isequal([0;1;0], round(final_pt))
                ends = ends + 2;
            end
            if isequal([0;0;1], round(final_pt))
                ends = ends + 3;
            end
        end
        result(array_length + 1 - i, j) = round(ends/double(samples_of_simplex));
    end
end

save("Results/PSD-1-antibiotic_KP-KD_fate", "result")