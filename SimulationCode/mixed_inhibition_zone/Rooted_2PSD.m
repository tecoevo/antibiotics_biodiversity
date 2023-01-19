% Here, we understand the effect of adding extra strains and extra
% antibiotics to the N=3, M=2 community which has a unique stable internal fixed point. 

addpath(./functions)

N=3;
M=2;

% metabolic costs of:
c_res = repelem(0.05,M);  % intrinsic resistance. uniformly chosen from [0.01, 0.06].
c_deg = 2.1*c_res;  % resistance via degradation
c_prod = 3*c_res; % production

K_D = 10;      %degrading area
K_P = 15;     %killing area

simplex_sampling_size = 500; %number of different initial conditions to supply to solver

P_o = [1, 0;
    0, 1;
    0, 0];
S_o = [0, 0;
    1, 0;
    0, 1];
D_o = [0, 1;
    0, 0;
    1, 0];
R_o = [0, 0;
    0, 0;
    0, 0];

cell_counter = 1;
communities_with_stable_coexistence = {};

% here we set the manner by which the community is made more complex

expand = 'A'; % determines in what direction the matrix is being expanded 
% - extra strain or extra antibiotic. can be set to 'S'(strain) or 'A'(antibiotic). 

expand_num = 1; % number of rows/columns by which the matrices is being expanded.

if expand == 'S'

    N = N + expand_num;
    g = repelem(1,N); % base fitness of each species

    load("NRcommunities/phenseq" + expand_num + "," + M + ".mat", 'phenotype_sequences')
    append_sequences = phenotype_sequences;

    for community_num = 1:size(append_sequences,2)

        sequence = append_sequences(:,community_num); % string of length N*M with entries in {1,2,3,4}
        P_add = zeros(expand_num,M);
        S_add = zeros(expand_num,M);
        D_add = zeros(expand_num,M);
        R_add = zeros(expand_num,M);

        for matrix_index = 1:expand_num*M %populating the matrices P, S, D, R
            % matrix_index will determine the matrix index (i,j)
            % where we record the phenotype of strain i with respect to antibiotic j

            phenotype = sequence(matrix_index,1);
            colnum = mod(matrix_index,M) + (M - mod(matrix_index,M))*isint(matrix_index/M);
            rownum = 1 + (matrix_index - modifiedmod(matrix_index,M))/M;
            if phenotype == 1
                P_add(rownum, colnum) = 1;
            elseif phenotype == 2
                S_add(rownum, colnum) = 1;
            elseif phenotype ==3
                D_add(rownum, colnum) = 1;
            else
                R_add(rownum, colnum) = 1;
            end
        end

        % if more strains are added, they need to be shuffled with respect
        % to antibiotics i.e., the columns of the additional matrices need to be
        % flipped in all possible combinations. For example, if the new
        % strain phenotype is PSD, the phenotype PDS is different. Writing
        % down the matrices and permuting the rows and columns shows that
        % these two additions do not lead to the same community.

        for flip_index = [0,1]
            if flip_index == 0
                P = [P_o; P_add]
                S = [S_o; S_add]
                D = [D_o; D_add]
                R = [R_o; R_add]
            elseif flip_index == 1
                % flipped append matrices
                P = [P_o; fliplr(P_add)]
                S = [S_o; fliplr(S_add)]
                D = [D_o; fliplr(D_add)]
                R = [R_o; fliplr(R_add)]
            end

            fixedpts = zeros(N,1);
            i=1; %to increase the size of the array fixedpts

            X = sym('X', [1 N], 'real');
            modeleqn = symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)-X;
            eqns = modeleqn == repelem(0,N);   % equations to solve

            for k = 1:simplex_sampling_size    % starting solver at different initial conditions
                solution = vpasolve(eqns, X,'Random', true);   % note: vpa might give multiple solutions
                fixedptarray = round(struct2array(solution) .* 1000)/1000; % each column is a fixed point
                for pt = fixedptarray
                    if is_internalfixedpt(pt) == 1 % record only "internal" fixed points i.e. those in the interior of the domain
                        fixedpts(:,i) = pt;
                        i=i+1;
                    end
                end
            end

            % Now we linearize about each fixed point and analyze stability.

            if isequal(fixedpts,repelem(0,N)')
                fprintf('Community number %d , %f: No internal fixed points. \n', community_num, flip_index);
            else
                disp(strcat(num2str(community_num),'_',num2str(flip_index)))
                fixedpts = unique(fixedpts', 'rows')' % neglect fixed points that have been repeated
            end

            Jacobian = jacobian(symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res), X);

            for ptcounter = 1:size(fixedpts,2) %iterating over fixed points of the system
                pt = fixedpts(:,ptcounter);
                if isequal(pt,repelem(0,N)')
                    continue
                end
                J_atspecificpoint = subs(Jacobian, X, pt');

                stability_parameter = double(vpa(log(max(abs(eig(J_atspecificpoint))))));

                if stability_parameter<0
                    communities_with_stable_coexistence{cell_counter}=strcat(num2str(community_num),'_',num2str(flip_index)); % record community number
                    cell_counter = cell_counter + 1;
                end
            end
        end
    end
    
elseif expand == 'A'

    M = M + expand_num;
    g = repelem(1,N); % base fitness of each species
    
    c_res = repelem(0.05,M);  % intrinsic resistance. uniformly chosen from [0.01, 0.06].
    c_deg = 2.1*c_res;  % resistance via degradation
    c_prod = 3*c_res;

    load("NRcommunities/phenseq" + N + "," + expand_num + ".mat", 'phenotype_sequences')
    append_sequences = phenotype_sequences;

    for community_num = [12, 14, 15]%1:size(append_sequences,2)
        community_num
        sequence = append_sequences(:,community_num); % string of length N*M with entries in {1,2,3,4}
        P_add = zeros(N,expand_num);
        S_add = zeros(N,expand_num);
        D_add = zeros(N,expand_num);
        R_add = zeros(N,expand_num);

        for matrix_index = 1:N*expand_num %populating the matrices P, S, D, R
            % matrix_index will determine the matrix index (i,j)
            % where we record the phenotype of strain i with respect to antibiotic j

            phenotype = sequence(matrix_index,1);
            colnum = mod(matrix_index,expand_num) + (expand_num - mod(matrix_index,expand_num))*isint(matrix_index/expand_num);
            rownum = 1 + (matrix_index - modifiedmod(matrix_index,expand_num))/expand_num;
            if phenotype == 1
                P_add(rownum, colnum) = 1;
            elseif phenotype == 2
                S_add(rownum, colnum) = 1;
            elseif phenotype ==3
                D_add(rownum, colnum) = 1;
            else
                R_add(rownum, colnum) = 1;
            end
        end
        
        % if more antibiotics are added, they need to be shuffled with respect
        % to the strains i.e., the rows of the additional matrices need to be
        % flipped in all possible combinations.
        for flip_index = [0,1,2,3,4,5]
            if flip_index == 0
                P = [P_o P_add]
                S = [S_o S_add]
                D = [D_o D_add]
                R = [R_o R_add]
            elseif flip_index == 1
                % flipped append matrices
                P = [P_o P_add([1 3 2],:)]
                S = [S_o S_add([1 3 2],:)]
                D = [D_o D_add([1 3 2],:)]
                R = [R_o R_add([1 3 2],:)]
            elseif flip_index == 2
                % flipped append matrices
                P = [P_o P_add([2 1 3],:)]
                S = [S_o S_add([2 1 3],:)]
                D = [D_o D_add([2 1 3],:)]
                R = [R_o R_add([2 1 3],:)]
            elseif flip_index == 3
                % flipped append matrices
                P = [P_o P_add([2 3 1],:)]
                S = [S_o S_add([2 3 1],:)]
                D = [D_o D_add([2 3 1],:)]
                R = [R_o R_add([2 3 1],:)]
            elseif flip_index == 4
                % flipped append matrices
                P = [P_o P_add([3 1 2],:)]
                S = [S_o S_add([3 1 2],:)]
                D = [D_o D_add([3 1 2],:)]
                R = [R_o R_add([3 1 2],:)]
            elseif flip_index == 5
                % flipped append matrices
                P = [P_o P_add([3 2 1],:)]
                S = [S_o S_add([3 2 1],:)]
                D = [D_o D_add([3 2 1],:)]
                R = [R_o R_add([3 2 1],:)]
            end


            fixedpts = zeros(N,1);
            i=1; %to increase the size of the array fixedpts

            X = sym('X', [1 N], 'real');
            modeleqn = symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)-X;
            eqns = modeleqn == repelem(0,N);   % equations to solve

            for k = 1:simplex_sampling_size    % starting solver at different initial conditions
                solution = vpasolve(eqns, X,'Random', true);   % note: vpa might give multiple solutions
                fixedptarray = round(struct2array(solution) .* 1000)/1000; % each column is a fixed point
                for pt = fixedptarray
                    if is_internalfixedpt(pt) == 1 % record only "internal" fixed points i.e. those in the interior of the domain
                        fixedpts(:,i) = pt;
                        i=i+1;
                    end
                end
            end

            % Now we linearize about each fixed point and analyze stability.

            if isequal(fixedpts,repelem(0,N)')
                fprintf('Community number %d , %f: No internal fixed points. \n', community_num, flip_index);
            else
                disp(strcat(num2str(community_num),'_',num2str(flip_index)))
                fixedpts = unique(fixedpts', 'rows')' % neglect fixed points that have been repeated
            end

            Jacobian = jacobian(symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res), X);

            for ptcounter = 1:size(fixedpts,2) %iterating over fixed points of the system
                pt = fixedpts(:,ptcounter);
                if isequal(pt,repelem(0,N)')
                    continue
                end
                J_atspecificpoint = subs(Jacobian, X, pt');

                stability_parameter = double(vpa(log(max(abs(eig(J_atspecificpoint))))));

                if stability_parameter<0
                    communities_with_stable_coexistence{cell_counter}=strcat(num2str(community_num),'_',num2str(flip_index)); % record community number
                    cell_counter = cell_counter + 1;
                end
            end
        end
    end
end

flag = false; 

if flag
    savename = "Results/Rooted/rooted2PSD_" + expand + "_" + expand_num + ".mat"
    save(savename,"communities_with_stable_coexistence")
end


celldisp(communities_with_stable_coexistence)


