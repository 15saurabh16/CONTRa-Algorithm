function [Final_trpl_list] = remove_redundant_tripoles(AllRootCents,AllLeaf1Cents,AllLeaf2Cents,AllDelVar,AllTsCorrMat,tau,true,chunk_size)
%% INPUT: location indices of the centroids of the three regions for all tripoles
% OUTPUT: location indice
PoleRedunMat = false(size(AllTsCorrMat));
PoleRedunMat(AllTsCorrMat>=tau) = true;

% addpath('~/expeditions/Saurabh/2016/');
% RedunDistMat = false(size(dist_mat));
% RedunDistMat(dist_mat<=beta) = true;
% 
% PoleRedunMat = PoleRedunMat & RedunDistMat;

Final_trpl_list = [];
[~,SortedDelVarInds] = sort(abs(AllDelVar),'descend');

AllChunkStInds = [1:chunk_size:length(SortedDelVarInds)];
AllChunkEndInds = AllChunkStInds+chunk_size-1;
for i = 1:length(AllChunkStInds)
%     i
    ChunkStSortedInd = AllChunkStInds(i);
    ChunkEndSortedInd = min(AllChunkEndInds(i),length(SortedDelVarInds));
    ChunkTrueInds = SortedDelVarInds(ChunkStSortedInd:ChunkEndSortedInd);
    
    ChunkRootCents = AllRootCents(ChunkTrueInds);
    ChunkLeaf1Cents = AllLeaf1Cents(ChunkTrueInds);
    ChunkLeaf2Cents = AllLeaf2Cents(ChunkTrueInds);
   
    ChunkDelVar = AllDelVar(ChunkTrueInds);
%     ChunkIndices = AllTrplIndices(:,ChunkTrueInds);
    
    if ~isempty(Final_trpl_list)
        SelectedRootCents = AllRootCents(Final_trpl_list);
        SelectedLeaf1Cents = AllLeaf1Cents(Final_trpl_list);
        SelectedLeaf2Cents = AllLeaf2Cents(Final_trpl_list);
%         SelectedStreng = AllStreng(Final_trpl_list);
%         SeletctedTrplIndices = AllTrplIndices(Final_trpl_list);

       [SelVsChnkRedunMat] = find_tripole_dist_mat(SelectedRootCents,SelectedLeaf1Cents,SelectedLeaf2Cents,ChunkRootCents,ChunkLeaf1Cents,ChunkLeaf2Cents,PoleRedunMat);
       RelChnkInds = find(~any(SelVsChnkRedunMat,1));
    else
       RelChnkInds = [1:length(ChunkTrueInds)];
    end
   
   
   if ~isempty(RelChnkInds)
       NewChnkTrueInds = ChunkTrueInds(RelChnkInds); 

       NewChunkRootCents = AllRootCents(NewChnkTrueInds);
       NewChunkLeaf1Cents = AllLeaf1Cents(NewChnkTrueInds);
       NewChunkLeaf2Cents = AllLeaf2Cents(NewChnkTrueInds);
       NewChunkDelVar = AllDelVar(NewChnkTrueInds);
%        NewChunkIndices = AllTrplIndices(:,NewChnkTrueInds);
       NewTrplRedunMat = find_tripole_dist_mat(NewChunkRootCents,NewChunkLeaf1Cents,NewChunkLeaf2Cents,NewChunkRootCents,NewChunkLeaf1Cents,NewChunkLeaf2Cents,PoleRedunMat);
       [NewTrplChnkList,~] = greedy_sel_nonredundant_tripoles(NewChunkDelVar,NewTrplRedunMat);
%        [NewTrplChnkList,~] = remove_redundant_tripoles_AND_dist_corr_no_trpl_index(NewChunkRootCents,NewChunkLeaf1Cents,NewChunkLeaf2Cents,NewChunkDelVar,PoleRedunMat,dist_mat,beta,tau,true);
       NewTrplList = NewChnkTrueInds(NewTrplChnkList);
       Final_trpl_list = [Final_trpl_list;NewTrplList];
   end
end



end