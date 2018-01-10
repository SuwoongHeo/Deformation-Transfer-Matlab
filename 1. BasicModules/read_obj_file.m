function [V,F,N] = read_vertices_and_faces_from_obj_file(filename)
  % Reads a .obj mesh file and outputs the vertex and face list
  % assumes a 3D triangle mesh and ignores everything but:
  % v x y z and f i j k lines
  % Input:
  %  filename  string of obj file's path
  %
  % Output:
  %  V  number of vertices x 3 array of vertex positions
  %  F  number of faces x 3 array of face indices
  %
  % Referred this site : http://www.alecjacobson.com/weblog/?p=917
  
  V = zeros(0,3);
  F = zeros(0,3);
  N = zeros(0,3);
  vertex_index = 1;
  normal_index = 1;
  face_index = 1;
  fid = fopen(filename,'rt');
  line = fgets(fid);
  while ischar(line)
    vertex = sscanf(line,'v %f %f %f');
    vertex_normal = sscanf(line, 'vn %f %f %f');
    face = sscanf(line,'f %d %d %d');
    face_long = sscanf(line,'f %d//%d %d//%d %d//%d',6);
    face_long_long = sscanf(line,'f %d/%d/%d %d/%d/%d %d/%d/%d',9);

    % see if line is vertex command if so add to vertices
    if(size(vertex)>0)
      V(vertex_index,:) = vertex;
      vertex_index = vertex_index+1;
   % see if line is simple face command if so add to faces
    elseif(size(face,1)==3)
      F(face_index,:) = face;
      face_index = face_index+1;
    % see if line is a face with normal indices command if so add to faces
    elseif(size(face_long,1)==6)
      % remove normal
      face_long = face_long(1:2:end);
      F(face_index,:) = face_long;
      face_index = face_index+1;
    % see if line is a face with normal and texture indices command if so add to faces
    elseif(size(face_long_long,1)==9)
      % remove normal and texture indices
      face_long_long = face_long_long(1:3:end);
      F(face_index,:) = face_long_long;
      face_index = face_index+1;
    elseif(size(vertex_normal,1)>0)
%       fprintf('Ignored: %s',line);
        N(normal_index, :) = vertex_normal;  
        normal_index = normal_index+1;
    else
        %TBD
    end

    line = fgets(fid);
  end
  fclose(fid);
end