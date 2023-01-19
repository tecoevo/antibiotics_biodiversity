% Here we examine the effect of metabolic cost on the stability of
% the fixed point of community number 79 from N=3, M=2.
% We assume that all species have an identical base fitness upon which
% corrections are made for production/degradation are made. These
% corrections are assumed to be in order c_prod > c_deg > c_res. 
% For convenience, we assume that c_prod = 3*c_res and c_deg = 2*c_res
% The non-redundant set of communities can be found in NRcommunities/phenseq3,2.mat

addpath("../functions/")

N = 3; % number of species
M = 2; % number of antibiotics

K_D = 5;      %degrading area
K_P = 10;     %killing area
g = repelem(1,N);

simplex_sampling_size = 200; %number of different initial conditions to supply to solver

load("NRcommunities/phenseq3,2.mat","phenotype_sequences")

% computing the matrices P,S,D,R

community_num = 79; 

sequence = phenotype_sequences(:,community_num); % string of length N*M with entries in {1,2,3,4}
P = zeros(N,M);
S = zeros(N,M);
D = zeros(N,M);
R = zeros(N,M);

for matrix_index = 1:N*M %populating the matrices P, S, D, R
    % matrix_index will determine the matrix index (i,j)
    % where we record the phenotype of strain i with respect to antibiotic j

    phenotype = sequence(matrix_index,1);
    colnum = mod(matrix_index,M) + (M - mod(matrix_index,M))*isint(matrix_index/M);
    rownum = 1 + (matrix_index - modifiedmod(matrix_index,M))/M;
    if phenotype == 1
        P(rownum, colnum) = 1;
    elseif phenotype == 2
        S(rownum, colnum) = 1;
    elseif phenotype ==3
        D(rownum, colnum) = 1;
    else
        R(rownum, colnum) = 1;
    end
end

maximum_cost_resistance = 0.18; % the PD phenotype takes the maximum hit of 0.9 i.e., 90 percent of its base rate
costs = 0.01:0.005:maximum_cost_resistance;

decay_rate = zeros(1,numel(costs));
X_1 = zeros(1,numel(costs));
X_2 = zeros(1,numel(costs));
X_3 = zeros(1,numel(costs));

for j = 1:numel(costs)
    try 
        
    c_res = costs(j);
    c_deg = 2.1*c_res;
    c_prod = 3*c_res;
    
    fixedpts = zeros(N,1);
    i=1; %to increase the size of the array fixedpts

    X = sym('X', [1 N], 'real');    
    modeleqn = symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)-X;
    eqns = modeleqn == repelem(0,N);   % equations to solve
    
    %solve for this value of metabolic costs
    for k = 1:simplex_sampling_size    % starting solver at different initial condition 
        solution = vpasolve(eqns, X,'Random', true);   % note: vpa might give multiple solutions
        fixedptarray = round(struct2array(solution) .* 1000)/1000; % each column is a fixed point
        for pt = fixedptarray
            if is_internalfixedpt(pt) == 1 % record only "internal" fixed points i.e. those in the interior of the domain
                fixedpts(:,i) = pt;
                i=i+1;
            end
        end
    end
    
    fixedpts = unique(fixedpts', 'rows')' % neglect fixed points that have been repeated
    
    X_1(j) = fixedpts(1,1);
    X_2(j) = fixedpts(2,1);
    X_3(j) = fixedpts(3,1);
    
     % stability of the fixed points
    Jacobian = jacobian(symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res),X);
    
    for ptcounter = 1:size(fixedpts,2) %iterating over fixed points of the system
        pt = fixedpts(:,ptcounter);
        J_atspecificpoint = subs(Jacobian, X, pt');
        decay_rate(j) = log(max(abs(eig(J_atspecificpoint))));
    end
    
    catch ME
        decay_rate(j) = NaN;
        X_1(j) = NaN;
        X_2(j) = NaN;
        X_3(j) = NaN;
    end
   
end

% plotting decay parameter with varying base growth rate
figure
plot(costs,decay_rate, "LineWidth",3)
title("Effect of metabolic costs on stability: stable if <0")
xlabel("Cost of resistance")
ylabel("Decay parameter")
xlim([0.005, 1.1*maximum_cost_resistance])
grid on
drawnow
saveas(gcf,'Figures/Metabolic_Cost/Metabolic_Cost-stability.png')


% plotting abundances of strains with varying cost of resistance
figure
plot(costs,X_1,'DisplayName','PD', "LineWidth",3)
hold on
plot(costs,X_2,'DisplayName','SP', "LineWidth",3)
hold on
plot(costs,X_3,'DisplayName','DS', "LineWidth",3)
hold on
title("Effect of metabolic cost on abundances at the fixed point")
xlabel("Cost of resistance")
ylabel("Abundances at fixed point")
xlim([0.005, 1.1*maximum_cost_resistance])
legend
grid on
drawnow
saveas(gcf,'Figures/Metabolic_Cost/Metabolic_Cost-abundance.png')



