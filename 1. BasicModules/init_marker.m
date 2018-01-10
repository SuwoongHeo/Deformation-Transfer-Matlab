function [ Marker ] = init_marker( VS, FS, VT, FT, FName )
%Initialize the user defiend marker
%   Graphical tool for selecting manual correspondence
%       Input
%           VS, FS : Source vertices and faces
%           VT, FT : Target vertices and faces
%           FName : String for saving the marker constraint
%       Output
%           Marker : Correspondence mapping

if exist(FName,'file')
    temp = open(FName);
    Marker = temp.Marker;
    return;    
end

prompt = 'Move cursor and press any key to continue(To stop[Y]/To remove[R]) : ';
error = 'Please grap the point and press any key to resume(To stop[Y]) : ';
%   Selection for source input
%% Marking source
figure; trimesh(FS, VS(:,1), VS(:,2), VS(:,3), ...
    'EdgeColor', 'none', 'FaceColor', [1 1 0], 'FaceLighting', 'phong'); 
light('Position',[0 0 1],'Style','infinite');
title('Marking in source mesh');
str = input('Align mesh and press any key to start', 's');
hold on; 
str = 'Any key';
cnt = 0;
while(1)    
    if strcmp(str, 'Y')
        break;
    elseif strcmp(str, 'R')        
        cnt = cnt-1;
        delete(handle1);
        str = input(prompt, 's');
    else
        obj1 = datacursormode();
        pt = getCursorInfo(obj1);
        if ~isempty(pt)
            cnt = cnt + 1;    
            ps = pt.Position;
            handle1 = scatter3(ps(1), ps(2), ps(3), 'r');
            S_ind(cnt, 1) = find( (VS(:,1) == ps(1)) & (VS(:,2) == ps(2)) & (VS(:,3) == ps(3)) );    
            str = input(prompt, 's');
        else
            str = input(error, 's');
        end
    end
end
hold off;

%% Marking Target
figure; trimesh(FT, VT(:,1), VT(:,2), VT(:,3), ...
    'EdgeColor', 'none', 'FaceColor', [1 1 0], 'FaceLighting', 'phong'); 
light('Position',[0 0 1],'Style','infinite');
title('Marking in target mesh');
hold on;
str = input('Align mesh and press any key to start', 's');
cnt = 0;
while(1)    
    if strcmp(str, 'Y')
        break;    
    elseif strcmp(str, 'R')        
        cnt = cnt-1;
        delete(handle2);
        str = input(prompt, 's');
    else
        obj2 = datacursormode();
        pt = getCursorInfo(obj2);
        if ~isempty(pt)
            cnt = cnt + 1;        
            ps = pt.Position;
            handle2 = scatter3(ps(1), ps(2), ps(3), 'g');
            T_ind(cnt, 1) = find( (VT(:,1) == ps(1)) & (VT(:,2) == ps(2)) & (VT(:,3) == ps(3)) );    
            str = input(prompt, 's');
        else
            str = input(error, 's');        
        end
    end
end
hold off;

Marker = [S_ind T_ind];
save(FName, 'Marker', 'S_ind', 'T_ind');
close all;

