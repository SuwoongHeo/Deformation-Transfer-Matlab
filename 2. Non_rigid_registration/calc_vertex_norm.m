function [ N ] = calc_vertex_norm(F, NF)
%Calculate vertex normal from adjacent face normal
%   F : Faces / NF : Normal of Faces
%   N : Normal of vertices
    len = max(max(F));
    N = zeros(len, 3);
    parfor i=1:len
        idx = find(F(:,1)==i | F(:,2) ==i | F(:,3) ==i);
        N(i,:) = sum(NF(idx,:))/length(idx);
        N(i,:) = N(i,:)/sqrt(sum(N(i,:).^2));
    end
end