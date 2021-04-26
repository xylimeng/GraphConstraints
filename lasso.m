function w = lasso(A, b, lambda, w)
% Lasso
p = size(A,1);
df = A * w - b;
eps = 1e-4;
finished = false(1, p);
while true
    for j = 1 : p
        wtmp = w(j);
        w(j) = soft(wtmp - df(j) / A(j,j), lambda / A(j,j));
        if w(j) ~= wtmp
            df = df + (w(j) - wtmp) * A(:, j); % update df
        end
    end
    index = (w == 0);
    finished(index) = (abs(df(index)) <= lambda);
    finished(~index) = (abs(df(~index) + lambda * sign(w(~index))) < eps);
    if finished
        break;
    end
end
end
%% Soft thresholding
function x = soft(x, lambda)
x = sign(x) * max(0, abs(x) - lambda);
end

