function [ corres ] = build_correspondence( VS, FS, VT, FT, maxind, thres , Name)
%Build correspondence using the proximity and face normals of source and
%target meshes
%   Input
%       VS : Defomred source mesh matched with target (nS x 3)
%       FS : Traingle indices of source mesh(mS x 3)
%       VT : Target mesh(nT x 3)
%       FT : Trinagle indices of target mesh(mT x 3)
%       maxind : Maximum correspondence
%       thres : Distance threshold for correspondence
%       Name : strings for speed up code
%   Output
%       corres : mT x # of correspondence for each triangles of target mesh

if nargin<5
    thres = 0.1;
    maxind = 10;
end
% %   Tempoerary code for speed up
if exist(Name, 'file')
    temp = open(Name);
    corres = temp.corres;    
    return;
end
fprintf('Building correspondene...\n');

[TS NS VS4 FS4]= v4_normal(VS, FS);
[TT NT VT4 FT4] = v4_normal(VT, FT);

VS_C = zeros(size(FS,1),3);
VT_C = zeros(size(FT,1),3);
%   Initial Triangle Correspondence
for i=1:length(FT)
    VT_C(i,:) = mean(VT(FT(i,:),:))';     %Centroids of target     
end
for i=1:length(FS)
    VS_C(i,:) = mean(VS(FS(i,:),:))';     %Centroids of source
end

S_tree = kd_tree(VS_C);
T_tree = kd_tree(VT_C);

corres1 = zeros(length(FT), maxind);
corres2 = zeros(length(FT), maxind);
templength = 0;
len = 0;
%% For source to target triangle coresspondence
% tic
fprintf('Source to Target correspondence..\n');
reverseStr=[];
rowlen = -1;    %##
for i=1:size(FS,1)  
    corresind = kd_query(T_tree, VS_C(i,:), thres, maxind);   
    corresind = corresind(corresind>0);
    len = length(corresind);
    corresind(sum(repmat(NS(i,:), [size(NT(corresind,:),1) 1]).*NT(corresind,:),2)>=pi/2) = 0;

    if ~isempty(corresind)
        for j = 1:len
            templength = max([rowlen+1 length(corres2(corres2(corresind(j),:)>0))]);
            rowlen = length(corres2(corres2(corresind(j),:)>0));
            corres2(corresind(j), rowlen+1) = i;
        end
    end
    
    if ~mod(i, 10000)
        msg = sprintf('Processed %d/%d', i, size(FS,1));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end
corres2 = corres2(:, 1:templength);

reverseStr = [];
fprintf('Target to Source correspondence..\n');

for i=1:size(FT,1)
    corresind = kd_query(S_tree, VT_C(i,:), thres, maxind);
    corresind = corresind(corresind>0);
    templength = max([len length(corresind)]);
    len = length(corresind);
    corresind(sum(repmat(NT(i,:), [size(NS(corresind,:),1) 1]).*NS(corresind,:),2)>=pi/2) = 0;
    corres1(i, 1:len) = corresind(1:len);

    if ~mod(i, 10000)
        msg = sprintf('Processed %d/%d', i, size(FT,1));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end
corres1 = corres1(:, 1:templength);
fprintf('\n');

tempcorres = [corres1 corres2];

corres = cell(size(FT,1),1);
for i=1:size(FT,1)
    temp = unique(tempcorres(i,:));
    corres{i} = temp(temp>0);
end

end


