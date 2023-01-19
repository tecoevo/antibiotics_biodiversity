function sequence = graph2seq(graph,N,M)

% converts the N-M bipartite graph to a sequence from which the matrices
% P,S,D,R can be computed. output is a column of length M*N

sequence = zeros(N*M,1);
i=1;

for strain = 1:N
    % look at phenotypes of this strain wrt all antibiotics
    for antibiotic = 1:M
        sequence(i,1) = graph.Edges.Weight(findedge(graph,strain, antibiotic+N));
        i=i+1;         
    end   
end  

end