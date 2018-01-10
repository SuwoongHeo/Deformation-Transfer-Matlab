function [ kd_tree ] = kd_tree( in )
%Implementation of KD - Tree 
%   Refer the slide which i made
%       Input : n x d Matrix (d is dimension)
%       Ouput : TBD

dim = size(in, 2);
ind = find(ones(size(in,1), 1));

kd_tree = partitioning(in, ind, dim);

end

function [ node ] = partitioning( in, ind, dim )
%Recursive method for tree creation
%   in : n x d Matrix for each partition
%   ind : index set of n input points
%   dim : dimension d
%   node : node of tree

    if size(in, 1) == 1
        node.pos = in;
        node.ind = ind;
        node.leaf = true;        
        return
    else
        node.leaf = false;
        subind = zeros(size(in,1), dim);
        for i=1:dim
            [temp subind(:,i)] = sort(in(:,i));
        end
    end
    clear temp;
    %max - min
    W = zeros(dim, 1);
    for i=1:dim
        W(i,1) = in(subind(end, i),i) - in(subind(1, i),i);
    end   
    
    [pline pdim] = max(W); % pdim is dimension of partioning
    
    node.axis = pdim; %line of partition
    node.pline =in(subind(1,pdim), pdim) + pline/2;
%     node.psubind = in(:,pdim)>node.pline;  %right = 1, left = 0
    psubind = in(:,pdim)>node.pline;  %right = 1, left = 0
    
%     l_in = in(~node.psubind, :);
%     l_ind = ind(~node.psubind);
    l_in = in(~psubind, :);
    l_ind = ind(~psubind);
%     for i=1:dim
%         [temp l_subind]= sort(l_in(:,i));
%     end
%     r_in = in(node.psubind, :);
%     r_ind = ind(node.psubind, :);
    r_in = in(psubind, :);
    r_ind = ind(psubind, :);
%     for i=1:dim
%         [temp r_subind(:,i)] = sort(r_in(:,i));
%     end    
    
    node.left = partitioning( l_in, l_ind, dim );
    node.right = partitioning( r_in, r_ind, dim );      
end

