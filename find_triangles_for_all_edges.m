function [All_triangles] = find_triangles_for_all_edges(All_edges,AllTsCorrMat,delta,RedunLocSets,trpl_type)

%% TRIPOLE SEARCH:
% Take care of both positive and negative tripoles 
% Remove kappa filtering for the third pole. 
% Output: triangles with three locations,strength,imp

All_triangles = [];
numLocs = size(AllTsCorrMat,1);

% METHOD1: 
TriplsForEdg = cell(length(All_edges),1);
for dipind = 1:size(All_edges,1) % PARFOR
    Tripoles_Dip = [];
    tripl_count = 0;
    
    loc1 = All_edges(dipind,1);Loc1Ts = AllTsCorrMat(loc1,:)';
    loc2 = All_edges(dipind,2);Loc2Ts = AllTsCorrMat(loc2,:)';
    Loc12_ts = zscore(Loc1Ts + Loc2Ts);
    edgecorr = AllTsCorrMat(loc1,loc2);   
    Loc1CorrMap  = AllTsCorrMat(:,loc1);
    Loc2CorrMap = AllTsCorrMat(:,loc2);
    
    
    ValidLocs = find(abs(Loc1CorrMap)<=abs(edgecorr) & abs(Loc2CorrMap)<=abs(edgecorr)); % edgecorr has to be strongest in magnitude since [loc1,loc2] is the super pair
    Loc1CorrsValid = Loc1CorrMap(ValidLocs);
    Loc2CorrsValid = Loc2CorrMap(ValidLocs);
    
    
    TrplStrengMapRoot3 = (Loc1CorrsValid + Loc2CorrsValid)./ sqrt(2*(1+edgecorr));%corr(Norm_ts',Loc12_ts);
    BestEdgMagn = max(abs(Loc1CorrsValid),abs(Loc2CorrsValid));
    DelVarMapRoot3 = (abs(TrplStrengMapRoot3).^2 - BestEdgMagn.^2);
    
    TrplStrengMapRoot1 = (Loc1CorrsValid + edgecorr)./ sqrt(2*(1+Loc2CorrsValid));%corr(Norm_ts',Loc12_ts);
    BestEdgMagn = max(abs(Loc1CorrsValid),abs(edgecorr));
    DelVarMapRoot1 = (abs(TrplStrengMapRoot1).^2 - BestEdgMagn.^2);
    
    TrplStrengMapRoot2 = (Loc2CorrsValid + edgecorr)./ sqrt(2*(1+Loc1CorrsValid));%corr(Norm_ts',Loc12_ts);
    BestEdgMagn = max(abs(Loc2CorrsValid),abs(edgecorr));
    DelVarMapRoot2 = (abs(TrplStrengMapRoot2).^2 - BestEdgMagn.^2);
    
    [DelVarMap,RootInd] = max([DelVarMapRoot1,DelVarMapRoot2,DelVarMapRoot3],[],2); % RootInd can take either 1,2, or 3 indicating which of the three poles is the root
    TrplStrengMap = zeros(size(DelVarMap));
    
    TrplStrengMap(RootInd==1) = TrplStrengMapRoot1(RootInd==1);
    TrplStrengMap(RootInd==2) = TrplStrengMapRoot2(RootInd==2);
    TrplStrengMap(RootInd==3) = TrplStrengMapRoot3(RootInd==3);
    
    BadValidInds1 =  find(DelVarMap<delta); % No Improvement
%     BadLocs2 = find(abs(Loc12CorrMap)<= abs(sigma)); % Poor strength are filtered out
    
    if strcmp(trpl_type,'pos')
       BadValidInds3 = find(TrplStrengMap<=0); % Since only positive tripoles are desired
       
    elseif strcmp(trpl_type,'neg')
       BadValidInds3 = find(TrplStrengMap>=0); % Since only negative tripoles are desired 
       
    else
       BadValidInds3 = []; % When both positive and negative tripoles are desired
    end
       
    BadValidInds = unique([BadValidInds1;BadValidInds3;loc1;loc2]); % Earlier BadLocs2 was also there
    GoodValidInds = setdiff([1:length(ValidLocs)]',BadValidInds);
    GoodRootInds = RootInd(GoodValidInds);
    GoodLocs = ValidLocs(GoodValidInds);
    GoodTrplStreng = TrplStrengMap(GoodValidInds);
    GoodDelVar = DelVarMap(GoodValidInds);
    
    numGoodLocs = length(GoodValidInds);
    Visited = logical(zeros(numLocs,1));
    
    % SANITY CHECK
    if min(GoodDelVar) < delta
        error(['Delta Threshold messing up for dipind:',num2str(dipind)]);
    end
    
    %% CODE 
    [sort_delvar,sort_goodvalidind] =sort(abs(GoodDelVar),'descend');

    sort_goodLoc = GoodLocs(sort_goodvalidind);
%     edge_inds = 1:length(Loc12CorrMap);
%     edge_ranks = zeros(length(Loc12CorrMap),1);
%     edge_ranks(sort_goodlocind) = edge_inds';
    loc_pointer = 0;
    num_seen = 0;   
    while true   
        loc_pointer = loc_pointer+1;
        if loc_pointer > numGoodLocs
            break;
        end
        if Visited(sort_goodLoc(loc_pointer))
            continue;
        end
        curr_goodvalidind = sort_goodvalidind(loc_pointer);
        loc3 = GoodLocs(curr_goodvalidind);% Alternatively = sort_goodLoc(loc_pointer); % AbsoluteLocind
%         Loc3_ts = zscore(mean(Norm_ts(:,loc3),2));
        trpl_streng = GoodTrplStreng(curr_goodvalidind);%corr(Loc12_ts,Loc3_ts);
        trpl_imp = GoodDelVar(curr_goodvalidind);
        edg1 = AllTsCorrMat(loc3,loc1);
        edg2 = AllTsCorrMat(loc3,loc2);
        edg3 = AllTsCorrMat(loc1,loc2);
        s_edge = max([abs(edg1),abs(edg2),abs(edg3)]);
%         trpl_imp = abs(trpl_streng/max(abs(edg1),abs(edg2)));
        Visited(RedunLocSets{loc3}) = true;
        ThisTrpl = [[loc2,loc3,loc1];[loc1,loc3,loc2];[loc1,loc2,loc3]]; 
        New_tripl = [ThisTrpl(GoodRootInds(curr_goodvalidind),:),trpl_streng,trpl_imp,s_edge,dipind];
        tripl_count = tripl_count+1;
        Tripoles_Dip(tripl_count,:) = New_tripl;                
    end
    
    TriplsForEdg{dipind} = Tripoles_Dip;
end

last_count = 0;
All_triangles = [];
for dipind= 1:size(All_edges,1)
    dipind
    sz_next_lot = size(TriplsForEdg{dipind},1);
    if sz_next_lot==0
        continue;
    end
    All_triangles(last_count+1:last_count+sz_next_lot,:) = TriplsForEdg{dipind}(:,1:5); % Total cols are 7, Excluding last two columns as they are only useful for debugging purposes
    last_count = last_count + sz_next_lot;
end
 
end