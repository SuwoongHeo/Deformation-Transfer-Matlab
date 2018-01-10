function [ M_P2 C_P2 ] = build_phase2( VS, FS, NS, VT, VTN, T_tree, marker, wc )
%Build phase 2 sparse matrix M_P2 ( Closest Valid Point term) with # of
%source vertices (nS), triangles (mS), target vertices(nT)
%   Input
%       VS     : Defomred Source mesh from previous step( nS x 3 )
%       FS     : Triangle index of source mesh (mS x 3)
%       NS     : Traingle normals of source mesh (mS x 3)
%       VT     : Target mesh ( nT x 3 )
%       VTN    : Vertex normals of target mesh ( nT x 3 )
%       T_tree : kd-tree of Target vertices
%       marker : marker constraint
%       wc     : Weight value
%   Output
%       M_P2   : (3 * nS) x (3 * (nS + mS)) big sparse matrix
%       C_P2   : (3 * nS) matrix

VSN = calc_vertex_norm(FS, NS);
S_size = size(VS,1);
valid_pt = zeros(S_size, 2);
C_P2 = zeros(3*S_size, 1);
reverseStr = [];
% tic
for j = 1:S_size
    if find(marker(:,1)==j)
        valid_pt(j, :) = [j marker(marker(:,1)==j,2)];
    else
%         valid_pt(j, :) = [j find_closest_validpt(T_tree, VS(j,:), VSN(j,:), VTN)];
        valid_pt(j, :) = [j find_closest_validpt(VS(j,:), VSN(j,:), VT, VTN)];
    end
    C_P2((1:3) + (j-1)*3, 1) = wc .* VT(valid_pt(j, 2),:)';
    if ~mod(j, 10000)
        msg = sprintf('Processed %d/%d', j, S_size);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end
fprintf('\n')
% toc

M_P2 = sparse((1:3*S_size)', (1:3*S_size)', repmat(wc, [3*S_size 1]), 3*S_size, 3*(length(VS)+length(FS)));

end


% function [ valid ] = find_closest_validpt(T_tree, spt, snormal, VTN)    
%     while(1)
%         [idx tpt] = kd_query(T_tree, spt);
%         if acos(snormal * VTN(idx,:)') < pi/2
%             valid = idx;
%             break;
%         else
%             T_tree = kd_delete(T_tree, tpt);
%         end
%     end       
% end

function [ valid ] = find_closest_validpt(spt, snormal, vpts, VTN)    
    d = sum((repmat(spt, [size(vpts,1), 1]) - vpts).^2,2);
    [~, ind] = sort(d);
    for i=1:length(d)        
        if acos(snormal * VTN(ind(i),:)') < pi/2
           valid = ind(i);
           break;
        end
    end
end
