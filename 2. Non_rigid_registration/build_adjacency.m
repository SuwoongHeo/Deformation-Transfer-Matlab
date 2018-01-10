function [ Adj_idx ] = build_adjacency( FS )
%Build up the Adjacency(3 Adjacency) matrix
%Note. 3 Adjacency means adjacent triangles which intersect with edges of
%reference triangle
%   Input
%       FS : Triangle indices of source mesh
%   Output
%       Adj_idx : # of triangle x 3 for 3-connectivity adjacency
%% 
Adj_idx = zeros(length(FS), 3); % Adjacent triangles for each edges
% reverseStr=[];
parfor i=1:length(FS)
    for j=1:3                
        idx = find(sum(FS==FS(i, j),2) & sum(FS==FS(i, mod(j,3)+1),2));
        if sum(idx~=i)
            Adj_idx(i,j) = idx(idx ~= i);
        end
    end
%     if ~mod(i, 10000)
%         msg = sprintf('Processed %d/%d', i, size(FS,1));
%         fprintf([reverseStr, msg]);
%         reverseStr = repmat(sprintf('\b'), 1, length(msg));
%     end
end
clear idx;

end

