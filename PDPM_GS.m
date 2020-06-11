%function [X, Y, TPR1, TNR1, E1, SP] = PDPM_GS(X0, Xr, Y0, S, p)
%
% output:
%
% X    :    reconstructed inverse covariance matrix by PDFPPA
% Y    :    intermediate variable
% obj_inner: inner loop objective function
% obj_outer: outer loop objective function
% X_Y:      the Frobenius norm of X-Y

function [X,Y,obj_inner,obj_outer,X_Y] = PDPM_GS(X0, Xr, Y0, S, p, Ty)
X = X0;
Y = Y0;
epsilon = 0;
%s = sum(sum(Xr~=0,2));
%mu = 1;
%gamma = 0.1;
%lambda = 0.3;  
%nu = 1;
mu = 1;
gamma = 0.01;
lambda = 0.02;  
nu = 1;

% mu = 0.1;
% gamma = 0.1 ;
% lambda = 0.001;
% nu = 10;
obj_inner = zeros(1,20);
obj_outer = zeros(1,20);
X_Y = zeros(1,20);
for l = 1:20
    l
    for k = 1:20
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
        Y1 = zeros(p);
        Y1 = (gamma/(gamma+nu))*X+(nu/(gamma+nu))*Y;
        count_on = 0;
        for i = 1:p
            for j = 1:p
                if Xr(i,j) == 0
                    Y1(i,j) = 0;
                end
            end
        end
        for i = 1:p
            count_on = count_on + Y1(i,i);
        end
        for i = 1:p
            Y1(i,i) = count_on/p;
        end
        count_off = zeros(6,1);
        n_off = zeros(6,1);
        for i = 1:(p-1)
            for j = (i+1):p
                if Xr(i,j) ~= 0
                    if Ty(i) == 1 && Ty(j) == 1
                        count_off(1,1) = count_off(1,1) + Y1(i,j);
                        n_off(1,1) = n_off(1,1) + 1;
                    end
                    if Ty(i) == 2 && Ty(j) == 2
                        count_off(2,1) = count_off(2,1) + Y1(i,j);
                        n_off(2,1) = n_off(2,1) + 1;
                    end
                    if Ty(i) == 3 && Ty(j) == 3
                        count_off(3,1) = count_off(3,1) + Y1(i,j);
                        n_off(3,1) = n_off(3,1) + 1;
                    end
                    if Ty(i) == 1 && Ty(j) == 2
                        count_off(4,1) = count_off(4,1) + Y1(i,j);
                        n_off(4,1) = n_off(4,1) + 1;
                    end
                    if Ty(i) == 2 && Ty(j) == 1
                        count_off(4,1) = count_off(4,1) + Y1(i,j);
                        n_off(4,1) = n_off(4,1) + 1;
                    end
                    if Ty(i) == 1 && Ty(j) == 3
                        count_off(5,1) = count_off(5,1) + Y1(i,j);
                        n_off(5,1) = n_off(5,1) + 1;
                    end
                    if Ty(i) == 3 && Ty(j) == 1
                        count_off(5,1) = count_off(5,1) + Y1(i,j);
                        n_off(5,1) = n_off(5,1) + 1;
                    end
                    if Ty(i) == 2 && Ty(j) == 3
                        count_off(6,1) = count_off(6,1) + Y1(i,j);
                        n_off(6,1) = n_off(6,1) + 1;
                    end
                    if Ty(i) == 3 && Ty(j) == 2
                        count_off(6,1) = count_off(6,1) + Y1(i,j);
                        n_off(6,1) = n_off(6,1) + 1;
                    end
                end
            end
        end
        for i = 1:(p-1)
            for j = (i+1):p
                if Xr(i,j) ~= 0
                    if Ty(i) == 1 && Ty(j) == 1
                         Y1(i,j) = count_off(1,1)/n_off(1,1);
                    end
                    if Ty(i) == 2 && Ty(j) == 2
                         Y1(i,j) = count_off(2,1)/n_off(2,1);
                    end
                    if Ty(i) == 3 && Ty(j) == 3
                         Y1(i,j) = count_off(3,1)/n_off(3,1); 
                    end
                    if Ty(i) == 1 && Ty(j) == 2
                        Y1(i,j) = count_off(4,1)/n_off(4,1);
                    end
                    if Ty(i) == 2 && Ty(j) == 1
                        Y1(i,j) = count_off(4,1)/n_off(4,1);
                    end
                    if Ty(i) == 1 && Ty(j) == 3
                        Y1(i,j) = count_off(5,1)/n_off(5,1);
                    end
                    if Ty(i) == 3 && Ty(j) == 1
                        Y1(i,j) = count_off(5,1)/n_off(5,1);
                    end
                    if Ty(i) == 2 && Ty(j) == 3
                        Y1(i,j) = count_off(6,1)/n_off(6,1);
                    end
                    if Ty(i) == 3 && Ty(j) == 2
                        Y1(i,j) = count_off(6,1)/n_off(6,1);
                    end
                end 
            end
        end
        for i = 2:p
            for j = 1:(i-1)
                Y1(i,j) = Y1(j,i);
            end
        end
        Y = Y1;        %the iterative scheme for updating Y
        %inner loop objective function
        obj_inner(1,k) = -log(det(X)) + trace(X) + norm(X,'fro') + norm(X-Y,'fro');
        %\|(X^{k+1},Y^{k+1})-(X^{k},Y^{k})\|_{F}^{2}
        %norm(X-X) + norm(Y-Y);
    end
    %outer loop objective function
    obj_outer(1,l) = -log(det(X)) + trace(X) + norm(X,'fro');
    %\|X^{k}-Y^{k}\|_{F}^{2}
    X_Y(1,l) = norm(X-Y,'fro');
    gamma = gamma * 2;   %increase gamma by *2 in each iteration
end






















 