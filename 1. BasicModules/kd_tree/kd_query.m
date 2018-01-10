function [ index pt ] = kd_query( kd_tree, in, ds, maxsize )
%Query the specified point
%   If the queried point exist, return that point
%   If not, return most adjacent node in kd-tree
%       input
%           kd_tree : kd-tree for query(struct)
%           in : 1 x d vector for query
%           ds : If specified, used for finding point within distance of ds
%           maxsize = max number of points of range search
%       output
%           index : index in kd_tree of querried output
%           pt : corresponding point 

if nargin<3
    ds = [];
    maxsize = 1;
end
p = [];
if isempty(ds)
    p = NNSearch( in, kd_tree , p, inf);
    index = p.ind;
    pt = p.pos;
else
    p.pos = inf*ones(maxsize, 4);
    p.ind = zeros(maxsize, 1);
    cnt = 1;
    p = RSearch( in, kd_tree, p, ds, cnt);
%     index = zeros(size(p));
%     index = zeros(maxsize,1);
%     pt = zeros(size(in,2)+1, size(p,2));
%     for i=1:size(p,2)
%         index(i) = p(i).ind;
%         pt(:,i) = p(i).pos';
%     end
%     p.pos = p.pos(p.pos(:,1)~=0,:);
%     p.ind = p.ind(p.ind(:,1)~=0);
    [temp tempind] = sort(p.pos(:,4));
    index = p.ind(tempind);
    pt = p.pos(tempind, :);
end
    


end

function [ p ] = NNSearch( q, n, p, w )
%Nearest Neighbor Search using kd-tree
%It is recursive steps
%   input
%       q : querried point
%       n : current node
%       p : current nearest node
%       w : distance threshold in NNS step
%   output
%       index : index of found node
%       pt : corresponding point 
if n.leaf
    if sqrt(sum((q - n.pos).^2)) < w
        p.ind = n.ind; p.pos = n.pos;
        return 
    else

    end
else
    %In first step,
    if w == inf
        if q(n.axis) <= n.pline;
            %Go to left
            p = NNSearch(q, n.left, p, w);
            w = sqrt(sum((p.pos - q).^2));
            %If pline is in side of q(n.axis) + w
            if q(n.axis) + w > n.pline
                p = NNSearch(q, n.right, p, w);
            end
        else
            %Go to right
            p = NNSearch(q, n.right, p, w);
            w = sqrt(sum((p.pos - q).^2));
            %If pline is in side of q(n.axis) - w
            if q(n.axis) - w <= n.pline
                p = NNSearch(q, n.left, p, w);
            end
        end                
    else
        %if w is finite
        if q(n.axis) - w <= n.pline
            p = NNSearch(q, n.left, p, w);
            w = sqrt(sum((p.pos - q).^2));
        end
        if q(n.axis) + w > n.pline
            p = NNSearch(q, n.right, p, w);
        end
    end 
end

end

function [ p cnt ] = RSearch( q, n, p, ds, cnt)
%Range Serach for kd-tree
%   Find the points within R ranges
%   Same as NN Search but we use ds rather than w
if n.leaf
   if cnt > length(p.pos)
       return;
   end
   dist = sqrt(sum((q - n.pos).^2));
   if dist < ds
       p.pos(cnt, :) = [n.pos dist];
       p.ind(cnt) = n.ind;
       cnt = cnt+1;
%        temp.ind = n.ind ; temp.pos = [n.pos dist];
%        p = [p temp];
   end
else
    %In first step,
    if q(n.axis) - ds <= n.pline
        [p cnt]= RSearch(q, n.left, p, ds, cnt);
%         w = sqrt(sum((p.pos - q).^2));
    end
    if q(n.axis) + ds > n.pline
        [p cnt]= RSearch(q, n.right, p, ds, cnt);
    end
end
end