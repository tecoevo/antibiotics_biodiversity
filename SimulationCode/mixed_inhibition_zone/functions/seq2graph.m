
function G = seq2graph(sequence,N,M)
% convert sequence in {1,2,3,4} to N-M bipartite graph
G = graph;
G = addnode(G,N+M);
G.Nodes.Type = [repelem({'S'},N) repelem({'A'},M)]';

for strain=1:N
    for antibiotic=1:M
        NewEdge = table([strain, N+antibiotic], sequence(antibiotic + (strain-1)*M, 1), 'VariableNames', {'EndNodes','Weight'});
        G = addedge(G,NewEdge);
    end    
end

end
