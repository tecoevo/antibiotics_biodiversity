N = 3; % number of species
M = 2; % number of antibiotics

K_D = 10;     %degrading area
K_P = 15;     %killing area

c_res = repelem(0.05, M);  % intrinsic resistance
c_deg = 2.1*c_res;  % resistance via degradation
c_prod = 3*c_res; % production

g = repelem(1,N);

% use the block below to generate the matrices P, S, D, R if you can identify the
% community only by its index from a phenseq file (found in
% ./NRcommunities)

% filename = "NRcommunities/phenseq" + N + "," + M + ".mat";
% load(filename, "phenotype_sequences")
%
% community_num = 122;
%
% sequence = phenotype_sequences(:,community_num); % string of length N*M with entries in {1,2,3,4}
%
% P = zeros(N,M);
% S = zeros(N,M);
% D = zeros(N,M);
% R = zeros(N,M);
%
% for matrix_index = 1:N*M %populating the matrices P, S, D, R
%     % matrix_index will determine the matrix index (i,j)
%     % where we record the phenotype of strain i with respect to antibiotic j

%     phenotype = sequence(matrix_index,1);
%     colnum = mod(matrix_index,M) + (M - mod(matrix_index,M))*isint(matrix_index/M);
%     rownum = 1 + (matrix_index - modifiedmod(matrix_index,M))/M;
%     if phenotype == 1
%         P(rownum, colnum) = 1;
%     elseif phenotype == 2
%         S(rownum, colnum) = 1;
%     elseif phenotype ==3
%         D(rownum, colnum) = 1;
%     else
%         R(rownum, colnum) = 1;
%     end
% end

P = [1 0 ; 0 1 ; 0 0 ];
S = [0 1 ; 1 0 ; 0 0 ];
D = [0 0 ; 0 0 ; 1 1 ];
R = [0 0 ; 0 0 ; 0 0 ];

% or just define your own community by explicit matrices P, S, D, R
% P = [1 0 0 0 0 1; 0 1 0 1 0 0; 0 0 1 0 1 0];
% S = [0 0 1 1 0 0; 1 0 0 0 1 0; 0 1 0 0 0 1];
% D = [0 1 0 0 1 0; 0 0 1 0 0 1; 1 0 0 1 0 0];
% R = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];

X = sym('X', [1 N]);
modeleqn = symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)-X;

% to draw phase portrait:
% set each of these to either 0 or 1 to choose which of these you want on
% the simplex plot
plot_speed = 0;  % greyscale shading of the velocity magnitudes
plot_arrow = 0;  % local direction of trajectory in a few randomly selected points 
plot_traj = 0;   % plot trajectories starting from random initial points

% to draw trajectory:
timesteps = 1000; % how long do you want to evolve the system for the trajectories above?

hold on
axis equal
grid on
% title("Phase Portrait")

% set axis limits
ylim([-0.2,1.2])
xlim([-0.2,1.2])

% plot vertices of the simplex
vertexa = project_onto_plane([1 0 0]);
scatter(vertexa(1), vertexa(2), 100, 'k', 'filled', 'HandleVisibility',"off")

vertexb = project_onto_plane([0 1 0]);
scatter(vertexb(1), vertexb(2), 100, 'k', 'filled', 'HandleVisibility',"off")

vertexc = project_onto_plane([0 0 1]);
scatter(vertexc(1), vertexc(2), 100, 'k', 'filled', 'HandleVisibility',"off")

% draw boundaries of the simplex by connecting vertices
plot([vertexa(1) vertexb(1)], [vertexa(2) vertexb(2)], 'k', "LineWidth",1, 'HandleVisibility',"off")
plot([vertexb(1) vertexc(1)], [vertexb(2) vertexc(2)], 'k', "LineWidth",1, 'HandleVisibility',"off")
plot([vertexc(1) vertexa(1)], [vertexc(2) vertexa(2)], 'k', "LineWidth",1, 'HandleVisibility',"off")

% calculate speed at each point

if plot_speed == 1
    coverage = 30000;
    simplex_mesh = UniformSampleSimplex(N,coverage);
    plane_mesh = zeros(2,coverage);
    j=1;
    for pt = simplex_mesh
        plane_mesh(:,j) = project_onto_plane(pt)';
        j = j+1;
    end

    speed_mesh = zeros(1, coverage);
    i=1;
    for pt = simplex_mesh
        speed_mesh(:,i) = norm(project_onto_plane(model(pt',N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)-pt'));
        i=i+1;
    end
    scatter(plane_mesh(1,:), plane_mesh(2,:), 10, speed_mesh, 'Filled', "HandleVisibility","off")
    colormap(brighten(repmat(linspace(1,0,coverage)', 1, 3),-0.9))
    colorbar('Ticks',[])
end

% plot fixed point for PS-SP-DD community
% uncomment the block below only if you are analyzing the PS-SP-DD community
% fixed = [0.368 0.368 0.264];
% fixed_2d = project_onto_plane(fixed);
% scatter(fixed_2d(1), fixed_2d(2), 'k', 'DisplayName', 'Unstable fixed point')

% plot arrows that show local direction of trajectory
if plot_arrow == 1
    arrow_coverage = 100;
    simplex_mesh = UniformSampleSimplex(N,arrow_coverage);
    arrow_mesh = zeros(2,arrow_coverage);
    j=1;
    for pt = simplex_mesh
        arrow_mesh(:,j) = project_onto_plane(pt)';
        j = j+1;
    end
    arrow_direction = zeros(2,arrow_coverage);
    i=1;
    for pt = simplex_mesh
        arrow_direction(:,i) = project_onto_plane(model(pt',N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)-pt');
        i=i+1;
    end
    quiver(arrow_mesh(1,:),arrow_mesh(2,:),arrow_direction(1,:),arrow_direction(2,:),0.5, "Color","g","DisplayName","Dynamics", "LineWidth",2)
end

% plot multiple trajectories

if plot_traj == 1
    %with legend
    initial_point = UniformSampleSimplex(N,1);
    initial_point2d = project_onto_plane(initial_point')'

    trajectory = TimeEvolution(initial_point, N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res, timesteps);

    trajectory2d = zeros(2,timesteps);
    for i = 1:size(trajectory,2)
        trajectory2d(:,i) = project_onto_plane(trajectory(:,i)')';
    end

    % plot initial point
    plot(initial_point2d(1),initial_point2d(2),'d',"MarkerSize",8,"Marker",'o',"MarkerFaceColor",'b', "DisplayName","Initial Conditions")

    % plot trajectory
    plot(trajectory2d(1,:),trajectory2d(2,:),"LineWidth",2,"Color","b", "DisplayName","Trajectories")

    for traj = 1:3 % how many extra trajectories do you want on the plot? Here, there will be 3 extra - so a total of 4 trajectories.
        % without legend
        initial_point = UniformSampleSimplex(N,1);
        initial_point2d = project_onto_plane(initial_point')';

        trajectory = TimeEvolution(initial_point, N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res, timesteps);

        trajectory2d = zeros(2,timesteps);
        for i = 1:size(trajectory,2)
            trajectory2d(:,i) = project_onto_plane(trajectory(:,i)')';
        end

        % plot initial point
        plot(initial_point2d(1),initial_point2d(2),'d',"MarkerSize",8,"Marker",'o',"MarkerFaceColor","b", "HandleVisibility","off")

        % plot trajectory
        plot(trajectory2d(1,:),trajectory2d(2,:),"LineWidth",2,"Color","b", "HandleVisibility","off")

    end
end


% labelling vertices

text(-0.1,0, "PS")
text(1.05, 0, "SP")
text(0.48, 0.93, "DD") 

if plot_traj == 1 || plot_arrow == 1
    legend
end

hold off
drawnow

function planept = project_onto_plane(simplexpt)
% projects a point on the 2-simplex onto the plane.
% assumes 3dpoint is a row of length 3.

planept = [(2*simplexpt(2) + simplexpt(3))/2, (0.5*sqrt(3))*simplexpt(3)];

end
