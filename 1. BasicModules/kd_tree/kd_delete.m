function [ kd_tree ] = kd_delete( kd_tree, in )
%Deletion operation of kd_tree
%   kd_tree : original kd_tree (Matlab structure type)
%   in      : 1 x n vector (n is dimension of vector)

kd_tree = deletion(kd_tree, in);

end

function [ n ] = deletion(n, in)

if n.leaf    
    n.leaf = -1;
else
    if in(n.axis) <= n.pline
        n.left = deletion(n.left, in);
        if n.left.leaf == -1
            n = rmfield(n, {'left', 'pline', 'axis'});
            n = n.right; 
            return;
        end
    end
    if in(n.axis) > n.pline
        n.right = deletion(n.right , in);
        if n.right.leaf == -1
            n = rmfield(n, {'right', 'pline', 'axis'});
            n = n.left;            
            return;            
        end
    end
end

end