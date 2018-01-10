function [ R, t, s, res ] = similarity_fitting( A, B )
% Fit the pointset on A to B using similarity transform
% First initial alignment step is follows
% O. Sorkine, "Least-Squares Rigid Motion Using SVD", 2016
% wi in this case, are set to be 1

if nargin ~= 2
    error('Missing parameters');
end

if size(A,1) > size(A,2)
    A = A';
end
if size(B,1) > size(B,2)
    B = B';
end

% Step1. Compute centroid for each data point
cent = [mean(A,2) mean(B,2)];

% Step2. center the vertices into origin
X = A - repmat(cent(:,1), [1 size(A,2)]);
Y = B - repmat(cent(:,1), [1 size(B,2)]);

% Step3. Compute covariance matrix, Since all component of W are 1,
% W=eye(size(A,2),size(A,2));
S = X*eye(size(A,2),size(A,2))*Y';

% Step4. Rotation vector is computed using SVD, see reference for detailed
% derivation. The final component of diagonal matrix is needed when rotation 
% is  totally reflection case
[U D V] = svd(S);
W = eye(size(V,1), size(V,1));
W(end,end) = det(V*U');
R = V*W*U';

% Step5. Compute optimal transform 
t = cent(:,2) - R*cent(:,1);

% s = 1.158;
s = 1.0; %*- Isotropic scale
% s = ones(3,1); %*- Anisotropic scale
if isreal(R)
    b0(1:4) = vrrotmat2vec(R);
    if ~isreal(b0)
        b0 = abs(b0);
    end
else
    disp(R);
    system('pause');
end
b0(5:7) = t';
b0(8) = s;
% b0(8:10) = s;

options = optimoptions('lsqnonlin','MaxFunctionEvaluations', 100000);
options.Algorithm = 'levenberg-marquardt';
[b, res, resi] = lsqnonlin(@(b) resSimXform( b,double(A),double(B) ),(double(b0)),[],[],options);
% [-Inf -Inf -Inf -Inf -Inf -Inf 0],[Inf Inf Inf Inf Inf Inf Inf]
n = size(A,2);
r = b(1:4);
t = b(5:7);
s = b(8);
% s = b(8:10);
R = vrrotvec2mat(r);
rot_A = diag(s) * R * A + repmat(t', 1, n);
res = sum(sqrt(sum((B - rot_A).^2)))/length(B);
end

% Compute residual of similarity X form between A, B
function result = resSimXform( b,A,B )

r = b(1:4);
t = b(5:7);
s = b(8);
% s = b(8:10);
n = size(A,2);

if ~isreal(r)
    a = 1;
end
R = vrrotvec2mat(r);
rot_A = diag(s) * R * A + repmat(t', 1, n);

% result = reshape((B - rot_A)', [1 prod(size(B))]);
result = sum(sum((B-rot_A).^2,2));

end

%% Najunil 
% % R * A - B
% H = B * (A)';
% [U,S,V] = svd(H);
% R = U*V';
% 
% if (det(R)<0)
%     fprintf('Reflection detected\n');
%     U(3,:)=-U(3,:);
%     R = U*V';
% end
% 
% D = diag(ones(3, 1));
% 
% centroid_A = mean(A');
% centroid_B = mean(B');
% t = -R * centroid_A' + centroid_B';
