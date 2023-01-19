function NRarray = NonRedundantCommunities(N,M)
% Generate all possible communities with N species and M antibiotics and 
% remove all the repeats - these are communities whose corresponding
% N-M bipartite graphs are isomorphic

phenotype_sequences = GenerateAllSequences(4,N*M);

NRgraphs = {};
i=1;

for sequence = phenotype_sequences
    querygraph = seq2graph(sequence,N,M);
    iso_checker = 1;
    query_edge_attr = querygraph.Edges.Weight;
    for ref = 1:numel(NRgraphs) 
        refgraph = NRgraphs{ref};
        ref_edge_attr = refgraph.Edges.Weight;
        % if the number of edges labelled 1, 2, 3 or 4 are not equal in the
        % 2 graphs, then they obviously cannot be isomorphic. The if-else
        % conditions below check these necessary (for an isomorphism to exist) 
        % conditions in an attempt to forego some isomorphism checks, which 
        % are expensive. 
        if sum(query_edge_attr==1) ~= sum(ref_edge_attr==1)
            iso_checker = iso_checker*1;
            continue
        elseif sum(query_edge_attr==2) ~= sum(ref_edge_attr==2)
            iso_checker = iso_checker*1;
            continue
        elseif sum(query_edge_attr==3) ~= sum(ref_edge_attr==3)
            iso_checker = iso_checker*1;
            continue
        elseif sum(query_edge_attr==4) ~= sum(ref_edge_attr==4)
            iso_checker = iso_checker*1;
            continue
        end
        % even if all of these quantities are equal, it is still possible that
        % the graphs are not isomorphic - equal edge numbers of each label is only a
        % necessary condition. Now we check isomorphism explicitly: 
        iso_compare = isisomorphic(querygraph, refgraph, 'EdgeVariables','Weight', 'NodeVariables', 'Type');
        if iso_compare == 0     % not isomorphic
            iso_checker = iso_checker*1;
        else                    % found an isomorphic graph
            iso_checker = 0;
            break
        end
    end
    % querygraph should be added to NRarray only if it is non-isomorphic 
    % to all graphs currently in NRarray i.e., if it is new. 
    % if there is at least 1 graph in NRarray isomorphic to querygraph, 
    % then we discard it since it is already represented. 
    if iso_checker == 1
        NRgraphs{i} = querygraph;
        i=i+1;
    end  
end

% convert set of nonredundant graphs to array of sequences in {1,2,3,4}. 
NRarray = zeros(N*M, numel(NRgraphs));

for graphid = 1:numel(NRgraphs) 
    graph = NRgraphs{graphid};
    NRarray(:,graphid) = graph2seq(graph,N,M);
end

end
