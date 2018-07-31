function [FinalTriangles] = CONTRa(AllTs,delta,kappa,tau,trpl_type)
%% INPUT:
% All_Ts: t x N matrix, where each column stores the time series of length t
% delta: A number in [0,1], Threshold on tripole jump
% kappa: A number in [0,1], Used in CONTRaFast to produce results faster
% with some compromize in completeness. Higher value leads to more computational time. 
% kappa = 1 is equivalent to CONTRaComplete. 

% tau:  [0,1]. A parameter that controls redundancy among tripoles in final
% output. 1 being most linent and considers every two tripoles to be
% different
% trpl_type: 
    %"pos" if only positive tripoles are desired, 
    % "neg" if only negative tripoles are desired, 
    %"all" if both pos and neg trpls are desired.
    
%% OUTPUT: 
% FinalTriangles: An m x 5 matrix, where m is the number of tripoles
% obtained. Each row indicates a tripole as [leaf1, leaf2, root, Signed Tripole Strength, Tripole Jump]

ALGOTIME = 0;

%% GENERATE CORRELATION MATRIX
AllTsCorrMat = corr(AllTs,AllTs);


%% FIND ZONES OF EXCLUSION
ExcludedZones =  find_redun_sets(AllTsCorrMat,kappa); % only used in CONTRaFast
RedunLocSets  = find_redun_sets(AllTsCorrMat,tau); % used to remove redundant tripoles

%% FIND CANDIDATE EDGES
tic;
[All_edges] = find_candidate_super_pairs(AllTsCorrMat,delta,ExcludedZones,trpl_type); 
disp(['All edges found:','Total No.',num2str(size(All_edges,1))]); 
DELT = toc;
ALGOTIME = ALGOTIME + DELT
 
NumSupPairs = length(All_edges)


%% FIND TRIPOLES FOR SUPER PAIR
AllGoodTr = find_triangles_for_all_edges(All_edges,AllTsCorrMat,delta,ExcludedZones,trpl_type);
DELT = toc;
ALGOTIME = ALGOTIME + DELT
NumTriangles = size(AllGoodTr,1)

%% SELECTING NON-REDUNDANT TRIPOLES
tic;
chunk_size = min(10000,size(AllGoodTr,1));
[Final_trpl_list] =  remove_redundant_tripoles(AllGoodTr(:,3),AllGoodTr(:,1),AllGoodTr(:,2),AllGoodTr(:,5),AllTsCorrMat,tau,true,chunk_size);
% [Final_trpl_list] = remove_redundant_tripoles_JAN_2017_trpl_index(AllGoodTr(:,3),AllGoodTr(:,1),AllGoodTr(:,2),AllGoodTr(:,5),AllGoodTrplIndices,DistMat,beta,tau,true,chunk_size);
DELT = toc;
FinalTriangles = AllGoodTr(Final_trpl_list,:);
ALGOTIME = ALGOTIME + DELT
end

