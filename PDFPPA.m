%function [X, Y, TPR1, TNR1, E1, SP] = PDFPPA(X0, Xr, Y0, S, p)
%
% output:
%
% X    :    reconstructed inverse covariance matrix by PDFPPA
% Y    :    intermediate variable
% E1   :    ELOSS value
% TPR1 :    TPR value
% TNR1 :    TNR value
% SP   :    number of non-zero elements of reconstructed inverse covariance matrix
%

function [X,Y] = PDFPPA(X0, Xr, Y0, S, p)
X = X0;
Y = Y0;
epsilon = 0.0001;
%epsilon = 0;
s = sum(sum(Xr~=0,2));
mu = 100;  %0.45
gamma = 100;
lambda = 0.001;
nu = 0.1;  %0.045
for l = 1:25
    for k = 1:15
       
        Z = mu*X+gamma*Y-S;
        [V,D] = eig(Z);
        for i = 1:p
            for j = 1:p
                T(i,j) = 0;
            end
        end
        for i = 1:p
            T(i,i)=(D(i,i)+sqrt(D(i,i)*D(i,i)+4*(lambda+gamma+mu)))/(2*(lambda+gamma+mu));
        end
        for i = 1:p
            T(i,i) = max(T(i,i),epsilon);
        end
        X = V*T*V';       %the iterative scheme for updating X
        Y1 = (gamma/(gamma+nu))*X+(nu/(gamma+nu))*Y;
        z1 = reshape(Y1,[],1);
        [z2, ind] = sort(abs(z1),'descend');
        ti = s + 1;
        z1(ind(ti:end)) = 0;
        Y1 = reshape(z1,p,p);
        Y = Y1;        %the iterative scheme for updating Y
    end
    gamma = gamma * 2;   %increase gamma by *2 in each iteration
end


SP = sum(sum(Y~=0,2)); %number of non-zero elements
[tmp,d] = eig(Y);
d = diag(d);
%compute objective value, ELOSS value
obj= sum(log(d)) - sum(sum(S*Y));
SY = inv(Xr)*Y; sx = eig(SY);
E1 = (sum(diag(SY)) - sum(log(sx)) - p)/p;
%Q1 = norm(SY - eye(p),'fro')/p;

%compute TPR value, TNR value
xr1 = 0; xr2 = 0;xr3 = 0;
for i = 1:p
    for j = 1:p
        if Xr(i,j) > 0
            xr1 = xr1+1;
        end
        if Xr(i,j) < 0
            xr2 = xr2+1;
        end
    end
end
xr3 = p*p - xr1 - xr2;
x_11 = 0; x_12 = 0; x_13 = 0;
for i = 1:p
    for j = 1:p
        if Y(i,j) > 0
            x_11 = x_11+1;
        end
        if Y(i,j) < 0
            x_12 = x_12+1;
        end
    end
end
x_13 = p*p - x_11 - x_12;
same3 = 0;
for i = 1:p
    for j = 1:p
        if Xr(i,j) ~= 0 && Y(i,j) ~= 0
            same3 = same3 + 1;
        end
    end
end
TPR1 = same3/(xr1+xr2);
same4 = 0;
for i = 1:p
    for j = 1:p
        if Xr(i,j) == 0 && Y(i,j) == 0
            same4 = same4 + 1;
        end
    end
end
TNR1 = same4/xr3;




 