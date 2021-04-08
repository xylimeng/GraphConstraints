function [Theta] = glasso(S,rho)
    [n,p] = size(S);
    max_iterations = 100;
    t = 1e-4;
    convergence_value = t * meanabs(S - diag(diag(S)));

    % initialise
    W_old = S + rho*eye(p);
    W = W_old;

    for round=1:max_iterations
        for j=p:-1:1
            i = j;

            W11 = W;
            W11(i,:) = [];  % remove ith row
            W11(:,j) = [];  % remove jth column
            w22 = W(i,j);

            s12 = S(:,j);
            s12(i,:) = [];
            
            A = W11^0.5;
            b = (W11^-0.5)*s12;
            
            beta = chenLasso(A,b,rho,1e2,1e-4);
            w12 = W11 * beta;

            W_left = W11(:,1:j-1);
            W_right = W11(:,j:p-1);
            W = [W_left w12 W_right];

            w12_row = [w12(1:j-1) ; w22 ; w12(j:p-1)]';
            W_above = W(1:i-1,:);
            W_below = W(i:p-1,:);
            W = [W_above ; w12_row ; W_below];
        end
        if meanabs(W - W_old) < convergence_value
           break; 
        end
        W_old = W;
    end
Theta = W^-1;



