function [T N V F] = v4_normal(Vert, Face)
%Convert 3 vertices representation to 4 vertices with representation of
%normal, 
%   input
%   Vert : n x 3 matrix
%   Face : m x 3 matrix
%   output
%   V : (m+n) x 3 matrix (vertices with added approximated normals)
%   T : n x 1 cell containing 3 x 3 matrices (representation of triangles)

f1 = Face(:,1); f2 = Face(:,2); f3 = Face(:,3);

e1 = Vert(f2,:) - Vert(f1,:);
e2 = Vert(f3,:) - Vert(f1,:);

c = cross(e1, e2, 2);
c_norm = sqrt(c(:,1).^2+c(:,2).^2+c(:,3).^2);
c_norm(c_norm==0)=1;

N = c./repmat(c_norm, [1 3]);
v4 = Vert(f1,:) + c./repmat(c_norm, [1 3]);

V = [Vert ;v4];
F = Face;
F(:, 4) = size(Vert,1)+find(F(:,3));
T = cell(size(F,1), 1);

for i = 1:size(F,1)
    T{i,1} = [(V(F(i,2), :) - V(F(i,1),:))' (V(F(i,3),:) - V(F(i,1),:))' (V(F(i,4),:) - V(F(i,1),:))'];
end

% quiver3(V(F(:,1),1), V(F(:,2),2), V(F(:,3),3), N(:,1), N(:,2), N(:,3));

end