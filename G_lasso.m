function [X, W] = G_lasso(S, lambda)
%% Graphical Lasso - Friedman et. al, Biostatistics, 2008
% Input:
%   S - sample covariance matrix
%   lambda - penalty parameter
% Output:
%   X - precision matrix    sigma^(-1)
%   W - covariance matrix   sigma
%%
p = size(S,1);  
W = S + lambda * eye(p);  %W=S+¦ËI
beta = zeros(p) - lambda * eye(p);   %¦Â=-¦ËI
eps = 1e-4;
finished = false(p);   
while true
    for j = 1 : p
        idx = 1 : p; idx(j) = [];
        beta(idx, j) = lasso(W(idx, idx), S(idx, j), lambda, beta(idx, j));
        W(idx, j) = W(idx,idx) * beta(idx, j);  %W=W*¦Â
        W(j, idx) = W(idx, j);
    end
    index = (beta == 0);
    finished(index) = (abs(W(index) - S(index)) <= lambda);
    finished(~index) = (abs(W(~index) -S(~index) + lambda * sign(beta(~index))) < eps);
    if finished
        break;
    end
end
X = zeros(p);
for j = 1 : p
    idx = 1 : p; idx(j) = [];
    X(j,j) = 1 / (W(j,j) - dot(W(idx,j), beta(idx,j)));
    X(idx, j) = -1 * X(j, j) * beta(idx,j);
end
% X = sparse(X);
end
