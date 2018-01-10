function [M C] = build_phase1(Adj_idx, E, FS4, VT4, ws, wi, marker)
% Create Sparse matrix M and vector C for phase 1 optimization.
%   Input arguments
%       Adj_idx : # of triangle x 3 matrix ( 3-connectivity )
%       E       : # of triangle cells (Each of 3 x 4 Elementary term
%       matrix)
%       FS4     : # of triangle x 4 matrix ( triangle index(4) matrix )
%       VT4     : # of vertices(target) x 3 matrix
%       ws, wi  : weight scalars
%       marker  : # of marker x 2 matrix
%   Output arguments
%       M       : Completed big sparse matrix( 9*(n_adj+length(FS)) x
%       3*(length(VT4) )
%       C       : Completed vector in rightside( 9*(n_adj+length(FS)) x 1 )
if nargin<4
    marker = zeros(1, 2);
end

n_adj = length(Adj_idx(:));
len_col = max(max(FS4));
I1 = zeros(9*n_adj*4, 3);   %index sets and value of U1
I2 = zeros(9*n_adj*4, 3);   %index sets and value of U2
I3 = zeros(9*length(FS4)*4, 3);
C1 = zeros(9*n_adj,1);
C2 = wi.*repmat(reshape(eye(3), [9 1]), [length(FS4) 1]);
%Sparse initialization 
%   Smothness term
% tic
reverseStr=[];
for i=1:length(FS4)                 %*- For ith triangle
    for j=1:3                       %*- For jth adjacent triangle
        if Adj_idx(i,j)  
            %Find the constraint vertices in triangle
            constid = zeros(2, 4);
            for k=1:3               %x, y, z
                %Find any marker constraint along ith triangle
                if sum(marker(:,1)==FS4(i,k))
                    constid(1, k) = k*(sum(marker(:,1)==FS4(i,k)));
                end
                %Find any marker constraint along jth adjcant triangle
                if sum(marker(:,1)==FS4(Adj_idx(i,j),k))
                    constid(2, k) = k*(sum(marker(:,1)==FS4(Adj_idx(i,j),k)));
                end
            end
            %Vectorized index of triangle U1(u1, u2, u3, u4) and U2(u1 ~)
            U1 = FS4(i,:); U2 = FS4(Adj_idx(i,j),:);
            for k=1:3   %x, y, z
                row = repmat((1:3)+(i-1)*3*3*3 + (j-1)*3*3 + (k-1)*3, [4 1]);
                col1 = repmat( (U1-1)*3 + k, [3 1])';
                val1 = ws.*E{i}';%*- value which corresponds to non-overlapping vertices
                if sum(constid(1,:))                    %Constrinat exist
                    C1((1:3)+(i-1)*3*3*3 + (j-1)*3*3 + (k-1)*3 , 1) = C1((1:3)+(i-1)*3*3*3 + (j-1)*3*3 + (k-1)*3 , 1) ...
                    -val1(constid(1,:)>0,:)' .* VT4(marker(marker==U1(constid(1,:)>0),2),k)';
                    val1(constid(1,:)>0,:) = 0;                    
                end
                col2 = repmat( (U2-1)*3 + k, [3 1])';
                val2 = -ws.*E{Adj_idx(i,j)}';%*- value which corresponds to overlapping vertices
                if sum(constid(2,:))                    %Constrinat exist
                    C1((1:3)+(i-1)*3*3*3 + (j-1)*3*3 + (k-1)*3 , 1) = C1((1:3)+(i-1)*3*3*3 + (j-1)*3*3 + (k-1)*3 , 1) ...
                    -val2(constid(2,:)>0,:)' .* VT4(marker(marker==U2(constid(2,:)>0),2),k)';
                    val2(constid(2,:)>0,:) = 0;
                end
                I1((1:12)+(i-1)*3*3*3*4 + (j-1)*3*3*4 + (k-1)*3*4,:) = [row(:) col1(:) val1(:)];
                I2((1:12)+(i-1)*3*3*3*4 + (j-1)*3*3*4 + (k-1)*3*4,:) = [row(:) col2(:) val2(:)];   
            end
        end
    end
    if i == 55507
        a = 1;
    end
%     if ~mod(i, 1)
%         msg = sprintf('Processed %d/%d', i, size(FS4,1));
%         fprintf([reverseStr, msg]);
%         reverseStr = repmat(sprintf('\b'), 1, length(msg));
%     end    
end
% fprintf('\n\n');
I1 = I1(I1(:,1)>0,:);
I2 = I2(I2(:,1)>0,:);

M1 = sparse(I1(:,1), I1(:,2), I1(:,3), 9*n_adj, 3*len_col);
M2 = sparse(I2(:,1), I2(:,2), I2(:,3), 9*n_adj, 3*len_col);
M3 = M1+M2;
clear I1 I2 col2 val2 constid ;
% toc

%   Identity term
%   Note. This term doesn`t need constraint (May be used with bigger wi
%   values)
% tic
for i=1:length(FS4)
    U1 = FS4(i,:);        
%     constid = zeros(1, 4);
%     for k=1:3
%         if sum(marker(:,1)==FS4(i,k))
%             constid(1, k) = k*(sum(marker(:,1)==FS4(i,k)));
%         end
%     end
    for k=1:3   %x, y, z
        row = repmat((1:3)+(i-1)*3*3 + (k-1)*3, [4 1]);
        col1 = repmat((U1-1)*3 + k, [3 1])';
        val1 = wi.*E{i}';     
%         if sum(constid(1,:))                    %Constrinat
%             C2((1:3)+(i-1)*3*3 + (k-1)*3 , 1) = C2((1:3)+(i-1)*3*3 + (k-1)*3 , 1) ...
%             -val1(constid(1,:)>0,:)' .* VT4(marker(marker==U1(constid(1,:)>0),2),k)';
%             val1(constid(1,:)>0,:) = 0;
%         end
        I3((1:12)+(i-1)*3*3*4 + (k-1)*3*4,:) = [row(:) col1(:) val1(:)];                
    end    
%     if ~mod(i, 10000)
%         msg = sprintf('Processed %d/%d', i, size(FS4,1));
%         fprintf([reverseStr, msg]);
%         reverseStr = repmat(sprintf('\b'), 1, length(msg));
%     end
end
% fprintf('\n\n');
clear row col1 val1;
% toc

M4 = sparse(I3(:,1), I3(:,2), I3(:,3));
C = [C1 ; C2];
M = [M3 ; M4];

clear M1 M2 M3 M4 I3;

end

