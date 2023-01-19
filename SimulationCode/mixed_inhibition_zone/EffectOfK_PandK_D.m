% Here we examine the stability behaviour of the fixed points of in the PD-SP-DS
% (N=3,M=2) community and the PDS-SPD-DS (N=3, M=3) community over the
% K_P - K_D parameter space.

% quantities describing the 2-antibiotic community will be henceforth
% indexed by "2" and quantities describing the 2-antibiotic community will be
% similarly indexed by "3". Common parameters will be left without an
% index.

N = 3; % number of species
M_2 = 2; % number of antibiotics
M_3 = 3;
g = repelem(1,N); % base fitness of each species

% metabolic costs of:
c_res_2 = repelem(1/20, M_2);  % intrinsic resistance. uniformly chosen from [0.01, 0.06].
c_deg_2 = 2.1*c_res_2;  % resistance via degradation
c_prod_2 = 3*c_res_2; % production

c_res_3 = repelem(1/20, M_3);  % intrinsic resistance. uniformly chosen from [0.01, 0.06].
c_deg_3 = 2.1*c_res_3;  % resistance via degradation
c_prod_3 = 3*c_res_3; % production

simplex_sampling_size = 200; %number of different initial conditions to supply to solver fmincon
timesteps = 600; %number of timesteps to evolve the system around an unstable fixed point

% setting the matrices P,S,D,R

P_2 = [1 0; 0 1; 0 0];
S_2 = [0 0; 1 0; 0 1];
D_2 = [0 1; 0 0; 1 0];
R_2 = [0 0; 0 0; 0 0];

P_3 = [1 0 0; 0 1 0; 0 0 1];
S_3 = [0 0 1; 1 0 0; 0 1 0];
D_3 = [0 1 0; 0 0 1; 1 0 0];
R_3 = [0 0 0; 0 0 0; 0 0 0];

% now consider different values of K_P and K_D

max_parameter_value = 2;

decay_rate_2 = zeros(max_parameter_value);

decay_rate_3 = zeros(max_parameter_value);

% we are looking at integer values of K_P and K_D each ranging from 1 to
% max_parameter_value. The integer value makes it convenient to index the stability matrix.
% For a particular K_P and K_D, stability_matrix(K_P, K_D) is 1 if the
% fixed point is stable and 0 if unstable.

% for the two-antibiotic community:

for K_P = 1:max_parameter_value
    K_P
    for K_D = 1:max_parameter_value
        fixedpts = zeros(N,1);
        i=1; %to increase the size of the array fixedpts

        X = sym('X', [1 N], 'real');
        modeleqn = symmodel(X,N,M_2,P_2,S_2,D_2,R_2,K_D,K_P,g,c_prod_2,c_deg_2,c_res_2)-X;
        eqns = modeleqn == repelem(0,N);   % equations to solve

        %solve for this value pair of K_P and K_D
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

        % stability of the fixed points
        try
            Jacobian = jacobian(symmodel(X,N,M_2,P_2,S_2,D_2,R_2,K_D,K_P,g,c_prod_2,c_deg_2,c_res_2),X);
            for ptcounter = 1:size(fixedpts,2) %iterating over fixed points of the system
                pt = fixedpts(:,ptcounter);
                J_atspecificpoint = subs(Jacobian, X, pt');
                decay_rate_2(max_parameter_value + 1 - K_P, K_D) = double(vpa(log(max(abs(eig(J_atspecificpoint))))));
            end
        catch ME % if Jacobian cannot be calculated i.e. internal fixed point doesn't exist
            decay_rate_2(max_parameter_value + 1 - K_P, K_D) = nan;
        end
    end
end

% for the three-antibiotic community:

for K_P = 1:max_parameter_value
    K_P
    for K_D = 1:max_parameter_value
        fixedpts = zeros(N,1);
        i=1; %to increase the size of the array fixedpts

        X = sym('X', [1 N], 'real');
        modeleqn = symmodel(X,N,M_3,P_3,S_3,D_3,R_3,K_D,K_P,g,c_prod_3,c_deg_3,c_res_3)-X;
        eqns = modeleqn == repelem(0,N);   % equations to solve

        %solve for this value pair of K_P and K_D
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

        % stability of the fixed points
        try
            Jacobian = jacobian(symmodel(X,N,M_3,P_3,S_3,D_3,R_3,K_D,K_P,g,c_prod_3,c_deg_3,c_res_3),X);
            for ptcounter = 1:size(fixedpts,2) %iterating over fixed points of the system
                pt = fixedpts(:,ptcounter);
                J_atspecificpoint = subs(Jacobian, X, pt');
                decay_rate_3(max_parameter_value + 1 - K_P, K_D) = double(vpa(log(max(abs(eig(J_atspecificpoint))))));
            end
        catch ME % if Jacobian cannot be calculated i.e. internal fixed point doesn't exist
            decay_rate_3(max_parameter_value + 1 - K_P, K_D) = nan;
        end
    end
end

mincolor = min([min(decay_rate_2), min(decay_rate_3)])
maxcolor = max([max(decay_rate_2), max(decay_rate_3)])
colormap_range = 1.2*max([abs(mincolor), abs(maxcolor)])

% saving raw data

save("Results/2PSD_vs_3PSD/KP_KD/2PSD_rawdata", "decay_rate_2")
save("Results/2PSD_vs_3PSD/KP_KD/3PSD_rawdata", "decay_rate_3")
 
% visualization: plotting effect on stability

f_2 = figure('visible','off')
h_2 = heatmap(num2cell(1:max_parameter_value), num2cell(max_parameter_value:-1:1), decay_rate_2)
h_2.GridVisible = "off";
h_2.CellLabelColor = "none";
h_2.MissingDataLabel = "Fixed point absent";
h_2.Colormap = colormap([autumn(64);flip(winter(64))]);
h_2.ColorLimits = [-colormap_range colormap_range];
title("Effect of K_P and K_D on dynamical behaviour near the fixed point")
xlabel("Production strength, K_P")
ylabel("Production strength, K_D")
saveas(f_2,'Results/2PSD_vs_3PSD/KP_KD/2PSD.png')

f_3 = figure('visible','off')
h_3 = heatmap(num2cell(1:max_parameter_value), num2cell(max_parameter_value:-1:1), decay_rate_3)
h_3.GridVisible = "off";
h_3.CellLabelColor = "none";
h_3.MissingDataLabel = "Fixed point absent";
h_3.Colormap = colormap([autumn(64);flip(winter(64))]);
h_3.ColorLimits = [-colormap_range colormap_range];
title("Effect of K_P and K_D on dynamical behaviour near the fixed point")
xlabel("Production strength, K_P")
ylabel("Production strength, K_D")
saveas(f_3,'Results/2PSD_vs_3PSD/KP_KD/3PSD.png')

delete(findall(0));
