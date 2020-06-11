%compute the interaction correlation between atoms from the estimated solution T 
%Input: 
%T: the estimated inverse covariance matrix T
%Xr: true inverse covariance matrix
%Ty: true atom type vector
%Output:
%rmse: RSME for one setting
%beta: the interaction correlation between atoms
%beta_true: the true interaction correlation between atoms
function [rmse,beta,beta_true] = post_procedure(T,p,Xr,Ty,K)
beta_true = zeros(6,1);
v = 1;
%T = Y;
count_beta = 0;
    for i = 1:p
        for j = (i+1):p
            if Xr(i,j) ~= 0
                if Ty(i,v) == 1 && Ty(j,v) == 1
                    beta_true(1,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
                if Ty(i,v) == 2 && Ty(j,v) == 2
                    beta_true(2,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
                if Ty(i,v) == 3 && Ty(j,v) == 3
                    beta_true(3,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
                if Ty(i,v) == 1 && Ty(j,v) == 2
                    beta_true(4,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
                if Ty(i,v) == 2 && Ty(j,v) == 1
                    beta_true(4,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
                if Ty(i,v) == 1 && Ty(j,v) == 3
                    beta_true(5,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
                if Ty(i,v) == 3 && Ty(j,v) == 1
                    beta_true(5,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
                if Ty(i,v) == 2 && Ty(j,v) == 3
                    beta_true(6,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
                if Ty(i,v) == 3 && Ty(j,v) == 2
                    beta_true(6,(v+1)/2) = Xr(i,j)*Xr(i,i);
                end
            end
        end
        %     for k = 1:6
        %         if beta(k,1) ~= 0
        %             count_beta = count_beta + 1;
        %         end
        %     end
        %     if count_beta == 6
        %         break;
        %     end
    end












beta = zeros(6,1);
v = 1;
%T = Y;
count_beta = 0;
    for i = 1:p
        for j = (i+1):p
            if Xr(i,j) ~= 0
                if Ty(i,v) == 1 && Ty(j,v) == 1
                    beta(1,(v+1)/2) = T(i,j)*T(i,i);
                end
                if Ty(i,v) == 2 && Ty(j,v) == 2
                    beta(2,(v+1)/2) = T(i,j)*T(i,i);
                end
                if Ty(i,v) == 3 && Ty(j,v) == 3
                    beta(3,(v+1)/2) = T(i,j)*T(i,i);
                end
                if Ty(i,v) == 1 && Ty(j,v) == 2
                    beta(4,(v+1)/2) = T(i,j)*T(i,i);
                end
                if Ty(i,v) == 2 && Ty(j,v) == 1
                    beta(4,(v+1)/2) = T(i,j)*T(i,i);
                end
                if Ty(i,v) == 1 && Ty(j,v) == 3
                    beta(5,(v+1)/2) = T(i,j)*T(i,i);
                end
                if Ty(i,v) == 3 && Ty(j,v) == 1
                    beta(5,(v+1)/2) = T(i,j)*T(i,i);
                end
                if Ty(i,v) == 2 && Ty(j,v) == 3
                    beta(6,(v+1)/2) = T(i,j)*T(i,i);
                end
                if Ty(i,v) == 3 && Ty(j,v) == 2
                    beta(6,(v+1)/2) = T(i,j)*T(i,i);
                end
            end
        end
        %     for k = 1:6
        %         if beta(k,1) ~= 0
        %             count_beta = count_beta + 1;
        %         end
        %     end
        %     if count_beta == 6
        %         break;
        %     end
    end
    
    
    
    
    
    
    
% %X = xlsread('beta_value.xlsx');
% X = beta;
% mean = zeros(6,1);
% for i = 1:6
%     for j = 1:K
%         mean(i,1) = mean(i,1) + X(i,j);
%     end
% end
% 
% for i = 1:6
%     mean(i,1) = mean(i,1)/K;
% end



%RMSE
%beta_true = [0.2,0.5,0.8, -0.3, -0.9, 0.7]';
rmse = zeros(6,1);
for i = 1:6
        rmse(i,1) = rmse(i,1) + (beta(i,1)- beta_true(i))^2;
end
for i = 1:6
    rmse(i,1) = sqrt(rmse(i,1));
end

