(* ::Package:: *)

% This code explores the space of all communities for a given N, M and saves the communities which have stable internal fixed points . The space of all communities can be % efficiently traversed by the algorithm outlined in Appendix \r ef{appendixA} of our manuscript . The results are saved as a cell array 'communities_with _stable _coexistence' which contains a community index k once for every fixed point this algorithm finds . unique(communities_with _stable _coexistence) will give a list of communities that have at least one stable internal fixed point .

addpath("../functions/")

N = 3; % number of species
M = 2; % number of antibiotics
disp("N and M set")
filename = "NRcommuniies/phenseq"+ N + "," + M + ".mat";
phenotype_sequences = NonRedundantCommunities(N,M);
disp('NR communities done')
save(filename, "phenotype_sequences")
disp("file saved")
exit
