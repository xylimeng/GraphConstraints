%Simulations in Section 4 of artical "Bayesian conditional autoregressive
%models with graph constraints for magnetic moment interaction in highentropy alloys"
clear all
close all
clc

p = 100;    %dimension
sparsity = 0.1;   %sparsity = 10p：p = 20 0.7; p = 50; 0.2; p = 100; 0.1; p = 200; 0.05
K = 50;   %number of different orders of atom type
mse = zeros(6,K); %RMSE
beta = zeros(6,K); %interaction correlation
beta_true = zeros(6,K); %true interaction correlation
count_beta  = zeros(6,1);
nonzero_true = zeros(1,K);
for a = 1:K
    a
    %geneate true inverse covariance matrix with different atom type vector
    [Xr,Xrt,Ty] = generate_true_revisit(p,sparsity,a);
    nonzero_true(1,a) = sum(sum(Xr~=0));
    samplenumber = 1; %sample number
    %generate different samples based on Xr
    [r_samples,S,X0,Y0] = generate_samples(a,p,Xrt,samplenumber);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % solve the problem by CARGO
    nu = p+1;
    c = 1;
    B_prior = c * eye(p);
    [B,T,obj_inner,obj_outer,X_Y]= CARGO(X0, Xr,Y0, S, p,Ty,nu,B_prior);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %solve the problem by Graphical Lasso
    % rho = 0.005;  %p = 50 0.01 p = 100 0.005 p = 20 0.001
    % [T,T1] = G_lasso(S, rho);
    % [T,T1] = graphical_lasso_1(S,rho,maxIt,tol)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %solve the problem by PDFPPA
    %[X,T] = PDFPPA(X0, Xr, Y0, S, p);
    [mse(:,a),beta(:,a),beta_true(:,a),count_beta_zero(:,a)] = post_procedure(T,p,Xr,Ty);
end
nonzero_mean = mean(nonzero_true);
%compute interaction correlations of atom magnetic moment over K settings
beta_aver = zeros(6,1);
beta_sum = zeros(6,1);
for i = 1:6
    for j = 1:K
        beta_sum(i,1) = beta_sum(i,1) + beta(i,j);
    end
end
for i = 1:6
    for j = 1:K
        count_beta(i,1) = count_beta(i,1) + count_beta_zero(i,j);
    end
end
for i = 1:6
    beta_aver(i,1) = beta_sum(i,1)/(K-count_beta(i,1));
end

%Compute mean square error 
mse_aver = zeros(6,1);
mse_sum = zeros(6,1);
for i = 1:K
    for j = 1:6
    mse_sum(j,1) = mse_sum(j,1) + mse(j,i);
    end
end
for i = 1:6
    mse_aver(i,1) = mse_sum(i,1)/(K-count_beta(i,1));
end
mse_mean = mean(mse_aver);
%standrad error of mse
std_mse = zeros(6,1);
for i = 1:6
    yy = mse(i,:);
    zz = yy(yy~=0);
    std_mse(i,1) = std(zz)/sqrt(K);
end
std_mean = mean(std_mse);

%plot Iterative step VS ||X-Y||_F and Iterative step VS Objective value: F(X)
%Implement only for CARGO
% x = linspace(1,35,35)
% %figure(1)
% plot(x,X_Y)
% xlabel('Iterative step','FontSize',15)
% ylabel('{||X-Y||_F}','FontSize',15);
% box off
% %hold on
% %figure(2)
% plot(x,obj_outer)
% xlabel('Iterative step','FontSize',15)
% ylabel('Objective value: F(X)','FontSize',15)
% box off


 
