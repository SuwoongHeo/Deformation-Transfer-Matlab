function [ VSP2 VT S_factor T_factor ] = non_rigid_registration( VS, FS, VT, FT, ws, wi, wc, marker, Name)
%Non-rigid registration of Deformation transfer
%   Input arguments
%      VS : Source vertices (n x 3 matrix)
%      FS : Triangle(Face) index of Source vertices (* x 3 matrix)
%      VT : Target vertices (m x 3 matrix)
%      FT : Triangle(Face) index of Target vertices (* x 3 matrix)
%      ws, wi, wc : weight paramters (scalar, scalar, nx1 vector)
%      marker : Constraint between Source and Target(Manual correspondence)
%      Name : strings for speed up code
%   Output arguements
%      VSP2 : Phase2 Optimized Result (n x 3 matrix)
%      VT   : VT with normalized coordinates (m x 3 matrix)
%      S_factor, T_factor : Normalizing factor for x,y,z coordinates(3 x 1)

if nargin <5
    ws = 1.0;
    wi = 0.001;
    wc = 1.0;    %1.0 to 5000.0 in paper..
    marker = [];
elseif nargin < 9
    Name = [];        
end
for i=1:length(marker)
    name{i} = i;
end
%Normalizing the coordinate of S and T to be insinde [-1,1]^3
S_factor = zeros(3,1);
T_factor = zeros(3,1);
S_size = size(VS, 1);
% T_size = size(VT, 1);

fprintf('Normalize Source & Target vert''s \n');
tmean = 0;
tstd = sqrt(2);
VS = normPts(VS, tmean, tstd);
VT = normPts(VT, tmean, tstd);

%   Tempoerary code for speed up
if exist(Name, 'file')
    fprintf('There already exist registred source file(DF_reg_phase2.mat)');
    fprintf('See through this file, if it is not the file you want, remove it');
    fprintf('and re-run this code\n');
    temp = open(Name);
    VSP2 = temp.VSP2;    
    return;
end

fprintf('Align Source to Target vert''s \n');
[R, t, s, res] = similarity_fitting(VT(marker(:,2),:), VS(marker(:,1),:));
VT = (VT*(s*R)' + repmat(t, length(VT), 1));

%   Converting to representation with,
%   Triangle(v2 - v1, v3 - v1, v4 - v1) :TS, TT
%   Normal(n1, n2, n3) : NS, NT
[TS NS VS4 FS4]= v4_normal(VS, FS);
[TT NT VT4 FT4] = v4_normal(VT, FT);


%%   Visualize
figure, fprintf('Visualize input meshes \n');
dispMesh(VS, FS, [.8 .8 .8], 0.8);hold on;
scatter3(VS(marker(:,1), 1), VS(marker(:,1), 2), VS(marker(:,1), 3), 'filled');
text(VS(marker(:,1), 1), VS(marker(:,1), 2), VS(marker(:,1), 3), name);
figure; 
dispMesh(VT, FT, [.8 0 .8], 0.8);hold on;
scatter3(VT(marker(:,2), 1), VT(marker(:,2), 2), VT(marker(:,2), 3), 'r', 'filled');
text(VT(marker(:,2), 1), VT(marker(:,2), 2), VT(marker(:,2), 3), name);
hold off;

fprintf('Building adjacency... \n');
Adj_idx = build_adjacency(FS);
fprintf('Solving Phase 1 optimization... \n');
E = build_elementary_cell(TS, length(FS));
[M C] = build_phase1(Adj_idx, E, FS4, VT4, ws, wi, marker);

fprintf('Apply phase 1 result \n');
%Construct Phase 1 result
VSP1 = M\C;
VSP1 = reshape(VSP1, [3 length(VSP1)/3])';
VSP1 = VSP1(1:S_size,:);

figure; title('Phase 1 registration results');
dispMesh(VSP1, FS, [.8 .8 .8], 0.8);
hold on;
scatter3(VSP1(marker(:,1), 1), VSP1(marker(:,1), 2), VSP1(marker(:,1), 3), 'filled');
text(VSP1(marker(:,1), 1), VSP1(marker(:,1), 2), VSP1(marker(:,1), 3), name);
hold on; dispMesh(VT, FT, [.8 0 .8], 0.8);
hold on;
scatter3(VT(marker(:,2), 1), VT(marker(:,2), 2), VT(marker(:,2), 3), 'r', 'filled');
text(VT(marker(:,2), 1), VT(marker(:,2), 2), VT(marker(:,2), 3), name);
hold off;


VTN = calc_vertex_norm(FT, NT);
T_tree = kd_tree(VT);
fprintf('Solving phase 2 optimization.. \n');
VSP2 = VSP1;
for i=1:length(wc)    
    ws = ws + (i-1)*wc(i)/100;
    [TS NS VS4 FS4]= v4_normal(VSP2, FS);
    E = build_elementary_cell(TS, length(FS));
    [M_P1 C_P1] = build_phase1(Adj_idx, E, FS4, VT4, ws, wi, marker);
    % mark.
    [M_P2 C_P2] = build_phase2(VSP2, FS, NS, VT, VTN, T_tree, marker, wc(i));
    
    M = [M_P1 ; M_P2];
    C = [C_P1 ; C_P2];
    
    VSP2 = M\C;
    VSP2 = reshape(VSP2, [3 length(VSP2)/3])';
    VSP2 = VSP2(1:S_size,:);
    
%     msg = sprintf('Processed %d/%d \n', i, length(wc));
%     fprintf(msg);    
end

figure; title('Phase 2 registration results');
trimesh(FS, VSP2(:, 1), VSP2(:, 2), VSP2(:, 3), ...
'EdgeColor', 'none', 'FaceColor', [1 1 0], 'FaceLighting', 'phong');
hold on;
trimesh(FT, VT(:, 1), VT(:, 2), VT(:, 3), ...
'EdgeColor', 'none', 'FaceColor', [0 1 1], 'FaceLighting', 'phong', 'facealpha', 0.4);
light('Position',[0 0 1],'Style','infinite');

figure; title('Phase 1 VS Phase 2');
trimesh(FS, VSP2(:, 1), VSP2(:, 2), VSP2(:, 3), ...
'EdgeColor', 'none', 'FaceColor', [1 1 0], 'FaceLighting', 'phong');
hold on;
trimesh(FS, VSP1(:, 1), VSP1(:, 2), VSP1(:, 3), ...
'EdgeColor', 'none', 'FaceColor', [0 1 1], 'FaceLighting', 'phong', 'facealpha', 0.4);
light('Position',[0 0 1],'Style','infinite');

save('DF_reg_phase2.mat', 'VSP2');

end


%%Mark
%     VSN = calc_vertex_norm(FS, NS);
%     
%     tic
%     for j = 1:S_size
%         if find(marker(:,1)==j)            
%             valid_pt(j, :) = [j marker(marker(:,1)==j,2)];
%         else
%             valid_pt(j, :) = [j find_closest_validpt(T_tree, VSP2(j,:), VSN(j,:), VTN)];            
%         end
%         C_P2((1:3) + (j-1)*3, 1) = wc(i) .* VT(valid_pt(j, 2),:)';   
%         if ~mod(j, 1000)
%             msg = sprintf('Processed %d/%d', j, S_size);
%             fprintf([reverseStr, msg]);
%             reverseStr = repmat(sprintf('\b'), 1, length(msg));
%         end
%     end
%     toc
%     fprintf('\n')
%     M_P2 = sparse((1:3*S_size)', (1:3*S_size)', repmat(wc(i), [3*S_size 1]), 3*S_size, 3*length(VS4));

%%
% function [ vs ] = opt_phase1( FS, TS_R, GLap, vs0, C )
% %Optimization function for non-rigid registrawtion
% %   자세한 설명 위치
% 
% vlength = max(max(FS));
% if nargin<4
%     vs0 = zeros(1, 3*vlength);
%     C = []
% end
% 
% options = optimset('TolFun', 1e-20, 'Display','iter');
% options.Algorithm = 'levenberg-marquardt';
% ub = 2*ones(1, 3*vlength);
% lb = -2*ones(1, 3*vlength);
% ub = [];
% lb = [];
% vs = lsqnonlin(@(vs) E_phase1( vs, FS, TS_R, GLap, 0.5, 0.005, C ), vs0, lb, ub, options);
% 
% end
% 
% function [ result ] = E_phase1( vs, FS, TS_R, GLap ,ws, wI, C )
% %
% %   TS_R is V mean traingles before deformation of source
% 
% %   Create triangle [v2-v1; v3-v1; v4-v1] after deformation
% vs = reshape(vs', [ 3 prod(size(vs'))/3 ])';
% [Td temp] = v4_normal(vs, FS);
% clear temp;
% T = cell(1,size(FS,1));
% % T = zeros(3,size(FS,1)*3);
% % T = zeros(9,size(FS,1));
% for i=1:size(FS,1)
%     T{1,i} = Td{i}/TS_R{i};
% %     T(1:3, ((i-1)*3+1):3*i) = Td{i}/TS_R{i};
% %     T(1:9, i) = reshape((Td{i}/TS_R{i})', [9 1]);
% end
% % 9 length vector 로 바꾸기
% 
% %   Smoothness term(E_S)
% E_S = 0;
% for i=1:size(FS,1)
%     ind = find(GLap(i,:) == -1);
%     for j=1:size(ind)
%         E_S = E_S + sum(sum((T{1,i} - T{1,ind(j)}).^2));
%     end
% end
% 
% %   Regularization term(E_I)
% E_I = 0;
% for i=1:size(FS,1)
%     E_I = E_I + sum(sum((T{1,i} - eye(3)).^2));
% end
% 
% result = ws*E_S + wI*E_I;
% if vs(prod(size(vs),2))
%     system('pause');
% end
% end

%% Adjacency code(Not used)
% SE = zeros(2, length(FS)+length(VS_norm)-2);    %Edge matrix
% tic
% cnt = 1;
% for i=1:length(FS)
%     for j=1:3
%         ind = mod(j, 3) + 1;
%         if ~sum(find(SE(2, SE(1,:)==FS(i, ind)...
%              )==FS(i,j)))
%             SE(: , cnt) = [FS(i, j) FS(i,ind)];
%             cnt = cnt+1;
%         end
%     end
% end
% 
% n_adj = length(SE);
% Adj_idx = zeros(2, n_adj); % Adjacent triangles for each edges
% for i=1:n_adj
%     [r1 c1] = ind2sub(size(FS), find(SE(1,i)==FS));
%     [r2 c2] = ind2sub(size(FS), find(SE(2,i)==FS));
%     corres = [];
%     for j=1:length(r1)
%         corres = [corres; r2(find(r1(j)==r2))];
%     end
%     Adj_idx(:,i) = corres;
% end
% Adj_idx(:,find(Adj_idx(1,:)==Adj_idx(2,:))) = 0;    %Eliminate boundary edges
% toc
% clear r1 c1 r2 c2;

%% Graph Laplacian code(Not used)

% fprintf('Construct Graph Laplacian.. \n');
% if exist('Face_S_GLap.mat')
%     temp = open('Face_S_GLap.mat');
%     S_GLap = temp.S_GLap;
%     clear temp;
% else
% S_incident = zeros(size(FS,1), size(FS,1), 'int8');
% S_GLap = spalloc(size(FS,1), size(FS,1), size(FS,1)*30);
% tic
% for i = 1:size(FS,1)    
%     for j = 1:size(VS_norm,2)
%             S_incident(i,((FS(:,1)==FS(i,j)) | (FS(:,2)==FS(i,j)) |...
%                 (FS(:,3)==FS(i,j)))) = -1;         
%     end       
%     ind = find(S_incident(i,:));
%     S_incident(i, i) = length(ind)-1;
%     S_GLap(i, ind) = double(S_incident(i, ind));
% end
% clear S_incident;
% toc
% end

% fprintf('Conducting phase 1 optimization.. \n');
% opt_phase1(FS, TS_R, S_GLap);