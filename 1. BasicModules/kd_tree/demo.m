close all;

A = [ 1 -1; -1 1; 1 1; -1 -1];

A_tree = kd_tree(A);
A_tree = kd_insert(A_tree, [0.5 0.5], 5);
test(A_tree);
A_tree = kd_delete(A_tree, [0.5 0.5]);

nnresult = kd_query(A_tree, [-0.8 -0.7]);

