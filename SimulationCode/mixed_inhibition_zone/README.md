This folder contains the scripts necessary for analysis of the mixed-inhibition zone model. This model is defined in the manuscript, which can be found in antibiotics-and-diversity/Manuscript/Working_draft/AGV_WorkingDraft.pdf

SCRIPTS:

The folder ``functions`` contains all the functions that are used in this analysis. 

get_NRcommunities.m can be used to get all non-redundant communities with N species and M antibiotics, which can be defined within the script.

SearchSpaceOfCommunities.m explores the space of communities found above and saves a cell array containing the indices (with multiplicity) in the .mat file of communities that contain stable fixed points.
These results are stored in Results/Search-Space/. The communities that these indices correspond to can be found by calling these indices in the appropriate (N,M) file from ./NRcommunities.

Rooted_2PSD.m and Rooted_3PSD.m contain analyses of the 2 base communities (see report) when strains evolve capabilities with respect to a new antibiotic, or when new strains are added.
These results are stored in Results/Extensions/.

PSD_KPKD_fate.m and PSD_KPKD_fate_analyse.m analyse the fate of the community that is a PSD motif on 1 antibiotic for a large range of K_P and K_D values. This is done to compare it to the spatial case, where there is nontrivial structure to the behaviour of this community over a similar parameter space. 
These results are stored in Results/PSD_compare_spatial/.

EffectOfK_PandK_D.m and EffectOfMetabolicCost.m examine the effects of these parameter on the community {PD, SP, DS}. See main text for an explanation of this notation. 

FOLDERS:

./functions contains the functions (written as individual .m files) necessary to run the above scripts.

./NRcommunities contains files that store non-redundant lists of communities for a few N,M value pairs.

./Results contains three subfolders - Extensions, SearchSpace, PSD_compare_spatial. 


