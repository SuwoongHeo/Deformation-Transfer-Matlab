function [ x, T] = normPts( x, m, s )
%Normalize coordinate values of x.
%   Note. Transformation is done in homogeneous coordinate
%   Note. k should be less than 3
%   x : n x k vector (n is number of elements, k is dimension in inhomogeneous coordinate)
%   mu : scalar which represents target center of normalization
%   s : scalar which represents average scale of normalization
n = size(x,1);
k = size(x,2);
T = eye(k+1);
mu = mean(x);
% mu = [mean(x(:,1)) mean(x(:,2)) mean(x(:,3))];
T(1:k, k+1) = (m-mu)';
% Average distance to mean location
mean_distance = mean(sum(sqrt((x - repmat(mu, [n 1])).^2),2));
% mean_distance = mean(sqrt((x(:,1)-mu(1)).^2+(x(:,2)-mu(2)).^2)+(x(:,3)-mu(3)).^2);
scale = s/mean_distance;
     
T = scale*T; T(k+1,k+1) = 1;
% T(1:3, 1:4) = sx*T(1:3, 1:4);
% Convert to homogeneous coordinate (Assuming w = 1)
x = [x ones(length(x),1)];
x = (T*x')';
x = x(:,1:k);    
end

