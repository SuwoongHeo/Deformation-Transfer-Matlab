function [ kd_tree ] = kd_insert( kd_tree, in , ind)
%Insertion operation of kd_tree
%   kd_tree : original kd_tree (Matlab structure type)
%   in      : 1 x n vector (n is dimension of vector)
%   Note. Should provide index of input point
%   Note. Before insert, you should add 'in' to original
%       matrix. 

kd_tree = insertion(kd_tree, in, ind);

end

function [ n ] = insertion(n, in, ind)

if n.leaf    
    %If leaf, part two values
    pt = [n.pos ; in];
    indset = [n.ind ind];
    n. leaf = false;
    dim = size(in,2);
    subind = zeros(2, dim);
    for i=1:dim
        [temp subind(:,i)] = sort(pt(:,i));
    end
    clear temp;
    W = zeros(size(in,2), 1);
    for i=1:dim
        W(i, 1) = abs(in(i) - n.pos(i));
    end
    n = rmfield(n, {'pos', 'ind'});
    [pline pdim] = max(W);
    
    n.axis = pdim;
    n.pline = pt(subind(1,pdim), pdim) + pline/2;
    psubind = pt(:,pdim)>n.pline;
    
    n.left.pos = pt(~psubind, :);
    n.left.ind = indset(~psubind);
    n.left.leaf = true;
    
    n.right.pos = pt(psubind, :);
    n.right.ind = indset(psubind);
    n.right.leaf = true;
else
    if in(n.axis) <= n.pline
        n.left = insertion(n.left, in, ind);
    end
    if in(n.axis) > n.pline
        n.right = insertion(n.right , in, ind);
    end
end

end

