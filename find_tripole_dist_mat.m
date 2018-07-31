function [TrplRedunMat,RedunRR,RedunL1L1,RedunL2L2,RedunL1L2,RedunL2L1] = find_tripole_dist_mat(Set1RootCents,Set1Leaf1Cents,Set1Leaf2Cents,Set2RootCents,Set2Leaf1Cents,Set2Leaf2Cents,Redun_mat)
n1 = length(Set1RootCents);
n2 = length(Set2RootCents);


RedunRR = Redun_mat(Set1RootCents,Set2RootCents);
RedunL1L1 = Redun_mat(Set1Leaf1Cents,Set2Leaf1Cents);
RedunL2L2 = Redun_mat(Set1Leaf2Cents,Set2Leaf2Cents);
RedunL1L2 = Redun_mat(Set1Leaf1Cents,Set2Leaf2Cents);
RedunL2L1 = Redun_mat(Set1Leaf2Cents,Set2Leaf1Cents);

RedunLeaves1= and(RedunL1L1,RedunL2L2);
RedunLeaves2 = and(RedunL1L2,RedunL2L1);

RedunLeaves = or(RedunLeaves1,RedunLeaves2);
TrplRedunMat = and(RedunRR,RedunLeaves);

end
