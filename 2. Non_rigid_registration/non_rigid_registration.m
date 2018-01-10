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
fprintf('Solving phase 2 optimization.. \n');
VSP2 = VSP1;
for i=1:length(wc)    
    ws = ws + (i-1)*wc(i)/100;
    [TS NS VS4 FS4]= v4_normal(VSP2, FS);
    E = build_elementary_cell(TS, length(FS));
    [M_P1 C_P1] = build_phase1(Adj_idx, E, FS4, VT4, ws, wi, marker);
    % mark.
    [M_P2 C_P2] = build_phase2(VSP2, FS, NS, VT, VTN, marker, wc(i));
    
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
