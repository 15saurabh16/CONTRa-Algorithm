 
function [All_edges] = find_candidate_super_pairs(Corr_mat,delta,ExcludedZones,trpl_type)
% Finds all the candidate super-pairs that can potentially form interesting tripoles
% for the given parameter setting


All_edges = [];
numLocs = size(Corr_mat,1);%D_srnn.latGrid*D_srnn.lonGrid;


% Calculate threshold on super pair absolute strength
All_s = [0:0.01:1];
Allresid = All_s.^2 + All_s.^3 -delta*(1-All_s);
s_thresh = All_s(find(Allresid>=0,1))% All_s(min(find(Allresid>=0)));
corr_vec = Corr_mat(:);
corr_vec(isnan(corr_vec)) = 0; 
[sortedg_val,sortedg_ind] = sort(abs(corr_vec),'descend');
sortedg_val_signed = corr_vec(sortedg_ind);
Visited = false(numLocs,numLocs);
num_edges = length(sortedg_val)
edge_pointer = 0;


while edge_pointer<num_edges
    edge_pointer = edge_pointer+1;

    if Visited(sortedg_ind(edge_pointer))
        continue;
    end   
    
    if sortedg_val(edge_pointer)<s_thresh 
        break;       
    end
    
    if (strcmp(trpl_type,'neg') && delta>0.0903 && sortedg_val_signed(edge_pointer)>-s_thresh)
        continue;
    end
    
    edge_streng = sortedg_val_signed(edge_pointer);
    edge_index = sortedg_ind(edge_pointer);
    [pt1, pt2] = ind2sub([numLocs,numLocs],edge_index);
    Visited(pt1,pt2) = true;
    Visited(pt2,pt1) = true;

    Nb1 = ExcludedZones{pt1};
    Nb2 = ExcludedZones{pt2};

    Visited(Nb1,Nb2) = true;
    Visited(Nb2,Nb1) = true;
    All_edges = [All_edges;[pt1,pt2,edge_streng]];
    % DEBUG
%     disp(['Edges so far = ',num2str(size(All_edges,1))]);
%     disp(['Last edge streng = ',num2str(edge_streng)]);
end



%% OLD CODE

% if strcmp(trpl_type,'neg') && delta > 0.0903
%     [sortedg_val,sortedg_ind] = sort(corr_vec,'ascend');
%     Visited = false(numLocs,numLocs);
%     num_good_pairs = length(find(corr_vec<s_thresh))
%     num_edges = length(sortedg_val);
%     edge_pointer = 0;
%     while edge_pointer<num_edges
%         edge_pointer = edge_pointer+1;
% 
%         if Visited(sortedg_ind(edge_pointer))
%             continue;
%         end
% 
%         if sortedg_val(edge_pointer)>s_thresh
%             break;       
%         end
%         edge_streng = sortedg_val(edge_pointer);
%         edge_index = sortedg_ind(edge_pointer);
%         [pt1, pt2] = ind2sub([numLocs,numLocs],edge_index);
%         Visited(pt1,pt2) = true;
%         Visited(pt2,pt1) = true;
% 
%         Nb1 = ExcludedZones{pt1};
%         Nb2 = ExcludedZones{pt2};
% 
%         Visited(Nb1,Nb2) = true;
%         Visited(Nb2,Nb1) = true;
%         All_edges = [All_edges;[pt1,pt2,edge_streng]];
%     end
%     
% else
%     
%     [sortedg_abs,sortedg_ind] = sort(abs(corr_vec),'descend'); % DIFF (descend,abs(),sortedg_abs)
%     sortedg_val = corr_vec(sortedg_ind);% EXTRA
%     Visited = false(numLocs,numLocs);
%     num_good_pairs = length(find(abs(corr_vec)>=s_thresh))% DIFF (abs() >=)
%     num_edges = length(sortedg_abs);
%     edge_pointer = 0;
%     while edge_pointer<num_edges
%         edge_pointer = edge_pointer+1;
% 
%         if Visited(sortedg_ind(edge_pointer))
%             continue;
%         end
% 
%         if sortedg_abs(edge_pointer)<s_thresh % DIFF (abs() <)
%             break;       
%         end
%         edge_streng = sortedg_val(edge_pointer);
%         edge_index = sortedg_ind(edge_pointer);
%         [pt1, pt2] = ind2sub([numLocs,numLocs],edge_index);
%         Visited(pt1,pt2) = true;
%         Visited(pt2,pt1) = true;
% 
%         Nb1 = ExcludedZones{pt1};
%         Nb2 = ExcludedZones{pt2};
% 
%         Visited(Nb1,Nb2) = true;
%         Visited(Nb2,Nb1) = true;
%         All_edges = [All_edges;[pt1,pt2,edge_streng]];
%     end
% 
% end



