function [ x nx ] = deformation_transfer( VS, FS, VT, FT, VS2, FS2, corres )
%Implementation of Deformation transfer
%   
%   자세한 설명 위치

lenFS = length(FS);
lenFT = length(FT);
SD = cell(lenFS, 1);

%%  Visualize
fprintf('Visualize input meshes \n');
set(gcf, 'Renderer', 'opengl');
figure; title('Source mesh');
trimesh(FS, VS(:, 1), VS(:, 2), VS(:, 3), ...
     'EdgeColor', 'none', 'FaceColor', [1 1 1], 'FaceLighting', 'phong');
light('Position',[0 0 1],'Style','infinite');
figure; title('Source deformed');
trimesh(FS2, VS2(:, 1), VS2(:, 2), VS2(:, 3), ...
     'EdgeColor', 'none', 'FaceColor', [1 1 1], 'FaceLighting', 'phong');
light('Position',[0 0 1],'Style','infinite');
figure; title('Target mesh');
trimesh(FT, VT(:, 1), VT(:, 2), VT(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', [0.5 1 0.5], 'FaceLighting', 'phong');
light('Position',[0 0 1],'Style','infinite');


%%
[TS NS VS4 FS4]= v4_normal(VS, FS);
[TS2 NS2 VS42 FS42] = v4_normal(VS2, FS2);
[TT NT VT4 FT4] = v4_normal(VT, FT);
clear NS N2 NT;
for i=1:lenFS
    SD{i} = TS2{i} / TS{i};
end

E = build_elementary_cell(TT, length(FT));

n_corres = sum(cellfun('length', corres));
n_non_corres = sum(cellfun('isempty',corres));
I = zeros(9*(n_corres+n_non_corres)*4, 3);
C = zeros(9*(n_corres+n_non_corres), 1);

% tic
fprintf('Transfer deformation..\n');
offset = 0;
offset2 = 0;
% reverseStr=[];
for i=1:lenFT
    lenCor = length(corres{i});
    Cor = corres{i};    
    U = FT4(i, :);
    if lenCor
        for j=1:lenCor            
            for k = 1:3 %x, y, z
                row = repmat((1:3)+ offset + (j-1)*3*3 + (k-1)*3, [4 1]);
                col1 = repmat( (U-1)*3 + k, [3 1] )';
                val1 = E{i}';
                I((1:12)+ offset2 + (j-1)*3*3*4 + (k-1)*3*4,:) = [row(:) col1(:) val1(:)];
            end        
            C((1:9)+ offset + (j-1)*9,1) = reshape(SD{Cor(j)}', [9 1]);
        end
        offset = offset + 3*3*lenCor;
        offset2 = offset2 + 3*3*lenCor*4;
    else
        for k=1:3
            row = repmat((1:3)+ offset + (k-1)*3, [4 1]);
            col1 = repmat( (U-1)*3 + k, [3 1] )';
            val1 = E{i}';
            I((1:12)+ offset2 + (k-1)*3*4,:) = [row(:) col1(:) val1(:)];
        end
        C((1:9)+ offset,1) = reshape(eye(3)', [9 1]);
        offset = offset + 3*3;
        offset2 = offset2 + 3*3*4;
    end
% 	if ~mod(i, 10000)
%         msg = sprintf('Processed %d/%d', i, size(FT,1));
%         fprintf([reverseStr, msg]);
%         reverseStr = repmat(sprintf('\b'), 1, length(msg));
%     end
end

% toc
% I = I(I(:,1)>0,:);

M = sparse(I(:,1), I(:,2), I(:,3), 9*(n_corres + n_non_corres), 3*length(VT4));
x = (M'*M)\(M'*C);
x = reshape(x, [3 length(x)/3])';
x = x(1:length(VT), :);
fprintf('Finsiehd\n');
[temp nx] = v4_normal(x, FT);
clear temp;

figure;
trimesh(FT, x(:, 1), x(:, 2), x(:, 3), ...
'EdgeColor', 'none', 'FaceColor', [0 1 1], 'FaceLighting', 'phong');
light('Position',[0 0 1],'Style','infinite');

end
%%  W.r.t. source triangle

% for i=1:lenFS
%     lenCor = length(corres{i});
%     Cor = corres{i};    
%     for j=1:lenCor
%         U = FT4(Cor(j), :);
%         for k = 1:3 %x, y, z
% %             row = repmat((1:3)+(i-1)*3*3*prelenCor + (j-1)*3*3 + (k-1)*3, [4 1]);
%             row = repmat((1:3)+ offset + (j-1)*3*3 + (k-1)*3, [4 1]);
%             col1 = repmat( (U-1)*3 + k, [3 1] )';
%             val1 = E{Cor(j)}';
%             I((1:12)+ offset2 + (j-1)*3*3*4 + (k-1)*3*4,:) = [row(:) col1(:) val1(:)];
%         end        
%         C((1:9)+ offset + (j-1)*9,1) = reshape(SD{i}', [9 1]);
%     end
%     offset = offset + 3*3*lenCor;
%     offset2 = offset2 + 3*3*lenCor*4;
% 	if ~mod(i, 10000)
%         msg = sprintf('Processed %d/%d', i, size(FS,1));
%         fprintf([reverseStr, msg]);
%         reverseStr = repmat(sprintf('\b'), 1, length(msg));
%     end
% end

