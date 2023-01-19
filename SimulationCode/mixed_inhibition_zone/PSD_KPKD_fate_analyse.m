N = 3; % number of species
M = 1; % number of antibiotics
g = repelem(1,N); % base fitness of each species
timesteps = 1000;

% metabolic costs of:
c_res = repelem(0.05, M);  %intrinsic resistance
c_prod = 3*c_res; %production
c_deg = 2.1*c_res;  %resistance via degradation

result = load("Results/PSD_compare_spatial/PSD-1-antibiotic_KP-KD_fate.mat");
result = cell2mat(struct2cell((result)));

step_size = 0.5;
max_parameter_value = length(result)*step_size;

AxisLabels = step_size:step_size:max_parameter_value;
% Convert each number in the array into a string
CustomAxisLabels = string(AxisLabels);
% Replace all but the fifth elements by spaces
CustomAxisLabels(mod(AxisLabels,5) ~= 0) = " ";

h = heatmap(result, 'Colormap', [0.172549, 0.482353, 0.713725; 0.843137, 0.29, 0.29; 0.670588, 0.85098, 1], 'ColorbarVisible', 'off');

% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomAxisLabels;
h.YDisplayLabels = flip(CustomAxisLabels);
h.Title = 'Fate of well-mixed PSD Populations';
h.XLabel = 'K_P';
h.YLabel = 'K_D';

saveas(gcf,'Results/PSD-1-antibiotic_KP-KD_fate.pdf')


