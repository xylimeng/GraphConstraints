%"Sparse inverse covariance matrix estimation with graph constraints"
%penalty decomposition proximal modification of Gauss-Seidle algorithm for 
%detecting atomic interactions of high entropy alloy 
clear all
close all
clc

%demo-simulations


%generate true (inverse) covariance matrix
p = 100;
sparsity = 0.14;  %set sparsity level
K = 1;   %numbelr of different orders of atom type 
rmse = zeros(6,K); %RMSE
beta = zeros(6,K); %interaction correlation
beta_true = zeros(6,K); %true interaction correlation
for a = 1:K 
%geneate true inverse covariance matrix with different atom type vector 
[Xr,Xrt,S1,Ty] = generate_true(p,sparsity,a);

%sample number
samplenumber = 1;
    
      
    %generate different samples based on Xr
    
    [r_samples,S,X0,Y0] = generate_samples(a,p,Xr,Xrt,samplenumber);    
    %Y0 = eye(p);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%solve the problem by penalty decomposition proximal modification of Gauss-Seidle algorithm

    [B,T,obj_inner,obj_outer,X_Y]= PDPM_GS(X0, Xr,Y0, S, p,Ty);
 
    %compute RMSE, beta
    [rmse(:,a),beta(:,a),beta_true(:,a)] = post_procedure(T,p,Xr,Ty,K);
    a = a + 1;
end

%compute interaction correlations of atom magnetic moment over K settings  
beta_aver = zeros(6,1);
beta_sum = zeros(6,1);
for i = 1:K
    beta_sum(1,1) = beta_sum(1,1) + beta(1,i);
    beta_sum(2,1) = beta_sum(2,1) + beta(2,i);
    beta_sum(3,1) = beta_sum(3,1) + beta(3,i);
    beta_sum(4,1) = beta_sum(4,1) + beta(4,i);
    beta_sum(5,1) = beta_sum(5,1) + beta(5,i);
    beta_sum(6,1) = beta_sum(6,1) + beta(6,i);
end

for i = 1:6
    beta_aver(i,1) = beta_sum(i,1)/K; 
end


rmse_aver = zeros(6,1);
rmse_sum = zeros(6,1);
for i = 1:K
    rmse_sum(1,1) = rmse_sum(1,1) + rmse(1,i);
    rmse_sum(2,1) = rmse_sum(2,1) + rmse(2,i);
    rmse_sum(3,1) = rmse_sum(3,1) + rmse(3,i);
    rmse_sum(4,1) = rmse_sum(4,1) + rmse(4,i);
    rmse_sum(5,1) = rmse_sum(5,1) + rmse(5,i);
    rmse_sum(6,1) = rmse_sum(6,1) + rmse(6,i);
end

for i = 1:6
    rmse_aver(i,1) = rmse_sum(i,1)/K; 
end

%plot
x=linspace(1,20,20)


% plot(x,X_Y)
% xlabel('iterative step')
% ylabel('the Frobenius norm of X-Y')

plot(x,obj_outer)
xlabel('iterative step')
ylabel('the value of objective function F(X)')
