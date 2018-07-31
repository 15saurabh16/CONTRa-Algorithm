function RedunLocSets = find_redun_sets(Corr_mat,thresh)
numLocs = size(Corr_mat,1);
RedunLocSets = cell(numLocs,1);

for pt = 1:numLocs
    RedunLocSets{pt} = find(Corr_mat(:,pt)>thresh);
end