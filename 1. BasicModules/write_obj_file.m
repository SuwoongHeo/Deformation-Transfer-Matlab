function write_obj_file( V, F, N, FName )
%Simple writing code of obj file
%   Input
%       V : Vertices
%       F : Faces
%       N : Normals
%       FName : Name of file
    fid = fopen(FName, 'wt');
    fprintf(fid, '#Deformation transfer result\n');
    fprintf(fid, '# This code is written by Suwoong Heo ; vereurer at gmail dot com\n');
    fprintf(fid, '# Original paper is wrriten by Robert W. Sumner, Jovan Popovic.\n');
    fprintf(fid, '# ''Deformation Transfer for Triangle Meshse (ACM ToG)''.\n');
    n_vertices = length(V);
    n_faces = length(F);
    fprintf(fid, '.obj format, %d vertices, %d triangles\n', n_vertices, n_faces);
    for i=1:n_vertices
        fprintf(fid, 'v   %f   %f   %f\n', V(i,1), V(i,2), V(i,3));
    end
    for i=1:n_vertices
        fprintf(fid, 'vn   %f   %f   %f\n', N(i,1), N(i,2), N(i,3));
    end
    for i=1:n_faces        
        fprintf(fid, 'f %d//%d %d//%d %d//%d\n', ...
            F(i,1), F(i,1), F(i,2), F(i,2), F(i,3), F(i,3));        
    end
end

