%compute the interaction correlation between atoms from the estimated solution T 
%Input: 
%T: the estimated inverse covariance matrix T
%Xr: true inverse covariance matrix
%Ty: true atom type vector
%Output:
%rmse: RSME for one setting
%beta: the interaction correlation between atoms
%beta_true: the true interaction correlation between atoms
function [mse,beta,beta_true,count_beta_zero] = post_procedure(T,p,Xr,Ty)
beta_true = zeros(6,1);
v = 1;
%T = Y;
count_true = zeros(6,1);
for i = 1:p
        for j = (i+1):p
            if Xr(i,j) ~= 0
                if Ty(i,v) == 1 && Ty(j,v) == 1
                    beta_true(1,1) = beta_true(1,1) - Xr(i,j)/Xr(i,i);
                    count_true(1,1) = count_true(1,1) + 1;
                end
                if Ty(i,v) == 2 && Ty(j,v) == 2
                    beta_true(2,1) = beta_true(2,1) - Xr(i,j)/Xr(i,i);
                    count_true(2,1) = count_true(2,1) + 1;
                end
                if Ty(i,v) == 3 && Ty(j,v) == 3
                    beta_true(3,1) = beta_true(3,1) - Xr(i,j)/Xr(i,i);
                    count_true(3,1) = count_true(3,1) + 1;
                end
                if Ty(i,v) == 1 && Ty(j,v) == 2
                    beta_true(4,1) = beta_true(4,1) - Xr(i,j)/Xr(i,i);
                    count_true(4,1) = count_true(4,1) + 1;
                end
                if Ty(i,v) == 2 && Ty(j,v) == 1
                    beta_true(4,1) = beta_true(4,1) - Xr(i,j)/Xr(i,i);
                    count_true(4,1) = count_true(4,1) + 1;
                end
                if Ty(i,v) == 1 && Ty(j,v) == 3
                    beta_true(5,1) = beta_true(5,1) - Xr(i,j)/Xr(i,i);
                    count_true(5,1) = count_true(5,1) + 1;
                end
                if Ty(i,v) == 3 && Ty(j,v) == 1
                    beta_true(5,1) = beta_true(5,1) - Xr(i,j)/Xr(i,i);
                    count_true(5,1) = count_true(5,1) + 1;
                end
                if Ty(i,v) == 2 && Ty(j,v) == 3
                    beta_true(6,1) = beta_true(6,1) - Xr(i,j)/Xr(i,i);
                    count_true(6,1) = count_true(6,1) + 1;
                end
                if Ty(i,v) == 3 && Ty(j,v) == 2
                    beta_true(6,1) = beta_true(6,1) - Xr(i,j)/Xr(i,i);
                    count_true(6,1) = count_true(6,1) + 1;
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

    for h = 1:6
            if count_true(h,1) == 0
                count_true(h,1) = 1;     %If some atom type doesn't exist.
            end
    end
    for j = 1:6
        beta_true(j,v) = beta_true(j,v)/ count_true(j,1);
%         beta_true(2,v) = beta_true(2,v)/ count_true(2,1);
%         beta_true(3,v) = beta_true(3,v)/ count_true(3,1);
%         beta_true(4,v) = beta_true(4,v)/ count_true(4,1);
%         beta_true(5,v) = beta_true(5,v)/ count_true(5,1);
%         beta_true(6,v) = beta_true(6,v)/ count_true(6,1);
    end
    
% T_orig = zeros(p);
% for i = 1:p
%     for j = 1:p
%         T_orig(i,j) = T(i,j)/T(i,i);
%     end
% end
% 
% 
% beta = zeros(6,1);
% v = 1;
% %T = Y;
% count_beta = zeros(6,1);
% count_beta_zero = zeros(6,1);   
%     for i = 1:p
%         for j = (i+1):p
%             if T(i,j) ~= 0
%                 if Ty(i,v) == 1 && Ty(j,v) == 1
%                     beta(1,1) = beta(1,1) + T_orig(i,j);
%                     count_beta(1,1) = count_beta(1,1) + 1;
%                 end
%                 if Ty(i,v) == 2 && Ty(j,v) == 2
%                     beta(2,1) = beta(2,1) + T_orig(i,j);
%                     count_beta(2,1) = count_beta(2,1) + 1;
%                 end
%                 if Ty(i,v) == 3 && Ty(j,v) == 3
%                     beta(3,1) = beta(3,1) + T_orig(i,j);
%                     count_beta(3,1) = count_beta(3,1) + 1;
%                 end
%                 if Ty(i,v) == 1 && Ty(j,v) == 2
%                     beta(4,1) = beta(4,1) + T_orig(i,j) ;
%                     count_beta(4,1) = count_beta(4,1) + 1;
%                 end
%                 if Ty(i,v) == 2 && Ty(j,v) == 1
%                     beta(4,1) = beta(4,1) + T_orig(i,j) ;
%                     count_beta(4,1) = count_beta(4,1) + 1;
%                 end
%                 if Ty(i,v) == 1 && Ty(j,v) == 3
%                     beta(5,1) = beta(5,1) + T_orig(i,j) ;
%                     count_beta(5,1) = count_beta(5,1) + 1;
%                 end
%                 if Ty(i,v) == 3 && Ty(j,v) == 1
%                     beta(5,1) = beta(5,1) + T_orig(i,j) ;
%                     count_beta(5,1) = count_beta(5,1) + 1;
%                 end
%                 if Ty(i,v) == 2 && Ty(j,v) == 3
%                     beta(6,1) = beta(6,1) + T_orig(i,j) ;
%                     count_beta(6,1) = count_beta(6,1) + 1;
%                 end
%                 if Ty(i,v) == 3 && Ty(j,v) == 2
%                     beta(6,1) = beta(6,1) + T_orig(i,j) ;
%                     count_beta(6,1) = count_beta(6,1) + 1;
%                 end
%             end
%         end
%     end
%     for h = 1:6
%             if count_beta(h,1) == 0
%                 count_beta(h,1) = 1;
%                 count_beta_zero(h,1) = count_beta_zero(h,1) + 1;
%             end
%     end




beta = zeros(6,1);
v = 1;
%T = Y;
count_beta = zeros(6,1);
count_beta_zero = zeros(6,1);   
    for i = 1:p
        for j = (i+1):p
            if T(i,j) ~= 0
                if Ty(i,v) == 1 && Ty(j,v) == 1
                    beta(1,1) = beta(1,1) - T(i,j)/T(i,i);
                    count_beta(1,1) = count_beta(1,1) + 1;
                end
                if Ty(i,v) == 2 && Ty(j,v) == 2
                    beta(2,1) = beta(2,1) - T(i,j)/T(i,i);
                    count_beta(2,1) = count_beta(2,1) + 1;
                end
                if Ty(i,v) == 3 && Ty(j,v) == 3
                    beta(3,1) = beta(3,1) - T(i,j)/T(i,i);
                    count_beta(3,1) = count_beta(3,1) + 1;
                end
                if Ty(i,v) == 1 && Ty(j,v) == 2
                    beta(4,1) = beta(4,1) - T(i,j)/T(i,i);
                    count_beta(4,1) = count_beta(4,1) + 1;
                end
                if Ty(i,v) == 2 && Ty(j,v) == 1
                    beta(4,1) = beta(4,1) - T(i,j)/T(i,i);
                    count_beta(4,1) = count_beta(4,1) + 1;
                end
                if Ty(i,v) == 1 && Ty(j,v) == 3
                    beta(5,1) = beta(5,1) - T(i,j)/T(i,i);
                    count_beta(5,1) = count_beta(5,1) + 1;
                end
                if Ty(i,v) == 3 && Ty(j,v) == 1
                    beta(5,1) = beta(5,1) - T(i,j)/T(i,i);
                    count_beta(5,1) = count_beta(5,1) + 1;
                end
                if Ty(i,v) == 2 && Ty(j,v) == 3
                    beta(6,1) = beta(6,1) - T(i,j)/T(i,i);
                    count_beta(6,1) = count_beta(6,1) + 1;
                end
                if Ty(i,v) == 3 && Ty(j,v) == 2
                    beta(6,1) = beta(6,1) - T(i,j)/T(i,i);
                    count_beta(6,1) = count_beta(6,1) + 1;
                end
            end
        end
    end
    for h = 1:6
            if count_beta(h,1) == 0
                count_beta(h,1) = 1;     %If some atom type doesn't exist.
                count_beta_zero(h,1) = count_beta_zero(h,1) + 1;
            end
    end
    for j = 1:6
        beta(j,v) = beta(j,v)/ count_beta(j,1);
%         beta(2,v) = beta(2,v)/ count_beta(2,1);
%         beta(3,v) = beta(3,v)/ count_beta(3,1);
%         beta(4,v) = beta(4,v)/ count_beta(4,1);
%         beta(5,v) = beta(5,v)/ count_beta(5,1);
%         beta(6,v) = beta(6,v)/ count_beta(6,1);
    end
mse = zeros(6,1);
for i = 1:6
    if beta(i,1) ~= 0
        mse(i,1) = mse(i,1) + (beta(i,1)- beta_true(i))^2;
    else
        mse(i,1) = 0;
    end
end
% for i = 1:6
%     rmse(i,1) = sqrt(rmse(i,1));
% end

