function dispMesh (shp, tl, color, alpha, normal, name, saveMesh)
% Note, shp has nx3 elements and tl has mx3 elements
% Note, If saveMesh activated, mesh is saved as obj file format
if nargin < 4
    alpha = 1;
end
if nargin < 5
    saveMesh = false;    
end

FV.vertices = shp;
FV.faces = tl;
if size(color,1) < 2
    patch(FV, 'facecolor', color, 'edgecolor', 'none', 'vertexnormalsmode', 'auto', 'FaceAlpha', alpha);
else
    patch(FV, 'FaceVertexCData',color,'facecolor', 'interp', 'edgecolor', 'none', 'vertexnormalsmode', 'auto', 'FaceAlpha', alpha);
end
camlight('headlight');
lighting phong;
material dull;
axis vis3d
axis equal;
% title(name);

if saveMesh
    write_obj_file(shp, tl, normal, sprintf('%s%s',name,'.obj'));
    frpintf('Mesh file %s is saved\n', name);
end



end

%% Deprecated (From Basel Face Model)
% function dispFace (shp, tl, rp)
% % 	if size(shp,1) > size(shp,2)
% %        shp = shp'; 
% %     end
%     if size(shp,2) == 1
%         shp = reshape(shp, [ 3 prod(size(shp))/3 ])'; 
%     end
% 
% 	set(gcf, 'Renderer', 'opengl');
% % 	fig_pos = get(gcf, 'Position');
% % 	fig_pos(3) = rp.width;
% % 	fig_pos(4) = rp.height;
% % 	set(gcf, 'Position', fig_pos);
% 	set(gcf, 'ResizeFcn', @resizeCallback);
% 
% 	mesh_h = trimesh(tl, shp(:, 1), shp(:, 3), shp(:, 2),...
%         'Edgecolor', 'none', 'FaceColor', 'interp', 'FaceVertexCData', ...
%         repmat(rp.color, [size(shp,1),1]),'FaceLighting', 'phong');
% 
% % 	set(gca, ...
% % 		'DataAspectRatio', [ 1 1 1 ], ...
% % 		'PlotBoxAspectRatio', [ 1 1 1 ], ...
% % 		'Units', 'pixels', ...
% % 		'GridLineStyle', 'none', ...
% % 		'Position', [ 0 0 fig_pos(3) fig_pos(4) ], ...
% % 		'Visible', 'off', 'box', 'off', ...
% % 		'Projection', 'perspective' ...
% % 		); 
% 	
% % 	set(gcf, 'Color', [ 0 0 0 ]); 
% 	view(180 + rp.phi * 180 / pi, 0);
% 
% 	material([.5, .5, .1 1  ])
% 	camlight('headlight');
% 	
% function resizeCallback (obj, eventdata)
% 	
% 	fig = gcbf;
% 	fig_pos = get(fig, 'Position');
% 
% 	axis = findobj(get(fig, 'Children'), 'Tag', 'Axis.Head');
% 	set(axis, 'Position', [ 0 0 fig_pos(3) fig_pos(4) ]);
	
