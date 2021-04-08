clear all
close all
clc
n_sample = 24;
n_training = 47;
n_insample = 24;
p = 24;
X_ = zeros(144,120);
X = zeros(p,120);
%dhcp
% X0 = csvread('magmom0.csv'); 
% X1 = csvread('magmom1.csv');
% X2 = csvread('magmom2.csv');
% X3 = csvread('magmom3.csv'); 
% X4 = csvread('magmom4.csv');
% X5 = csvread('magmom5.csv');
% X6 = csvread('magmom6.csv');
% X7 = csvread('magmom7.csv');
% X8 = csvread('magmom8.csv');
% X9 = csvread('magmom9.csv');
% X10 = csvread('magmom10.csv');
% X11 = csvread('magmom11.csv');
% X12 = csvread('magmom12.csv');
% X13 = csvread('magmom13.csv');
% X14 = csvread('magmom14.csv');
% X15 = csvread('magmom15.csv');
% X16 = csvread('magmom16.csv');
% X17 = csvread('magmom17.csv');
% X18 = csvread('magmom18.csv');
% X19 = csvread('magmom19.csv');
% X20 = csvread('magmom20.csv');
% X21 = csvread('magmom21.csv');
% X22 = csvread('magmom22.csv');
% X23 = csvread('magmom23.csv');


%fcc
% X0 = csvread('magmom0-fcc.csv'); 
% X1 = csvread('magmom1-fcc.csv');
% X2 = csvread('magmom2-fcc.csv');
% X3 = csvread('magmom3-fcc.csv'); 
% X4 = csvread('magmom4-fcc.csv');
% X5 = csvread('magmom5-fcc.csv');
% X6 = csvread('magmom6-fcc.csv');
% X7 = csvread('magmom7-fcc.csv');
% X8 = csvread('magmom8-fcc.csv');
% X9 = csvread('magmom9-fcc.csv');
% X10 = csvread('magmom10-fcc.csv');
% X11 = csvread('magmom11-fcc.csv');
% X12 = csvread('magmom12-fcc.csv');
% X13 = csvread('magmom13-fcc.csv');
% X14 = csvread('magmom14-fcc.csv');
% X15 = csvread('magmom15-fcc.csv');
% X16 = csvread('magmom16-fcc.csv');
% X17 = csvread('magmom17-fcc.csv');
% X18 = csvread('magmom18-fcc.csv');
% X19 = csvread('magmom19-fcc.csv');
% X20 = csvread('magmom20-fcc.csv');
% X21 = csvread('magmom21-fcc.csv');
% X22 = csvread('magmom22-fcc.csv');
% X23 = csvread('magmom23-fcc.csv');



%hcp
X0 = csvread('magmom0-hcp.csv'); 
X1 = csvread('magmom1-hcp.csv');
X2 = csvread('magmom2-hcp.csv');
X3 = csvread('magmom3-hcp.csv'); 
X4 = csvread('magmom4-hcp.csv');
X5 = csvread('magmom5-hcp.csv');
X6 = csvread('magmom6-hcp.csv');
X7 = csvread('magmom7-hcp.csv');
X8 = csvread('magmom8-hcp.csv');
X9 = csvread('magmom9-hcp.csv');
X10 = csvread('magmom10-hcp.csv');
X11 = csvread('magmom11-hcp.csv');
X12 = csvread('magmom12-hcp.csv');
X13 = csvread('magmom13-hcp.csv');
X14 = csvread('magmom14-hcp.csv');
X15 = csvread('magmom15-hcp.csv');
X16 = csvread('magmom16-hcp.csv');
X17 = csvread('magmom17-hcp.csv');
X18 = csvread('magmom18-hcp.csv');
X19 = csvread('magmom19-hcp.csv');
X20 = csvread('magmom20-hcp.csv');
X21 = csvread('magmom21-hcp.csv');
X22 = csvread('magmom22-hcp.csv');
X23 = csvread('magmom23-hcp.csv');


% Ni 1 Co 2 Cr 3 Fe 4
tic
for j=1:2
    X(:,j) = X0(:,j);
end
for j=3:4
    X(:,j) = X1(:,j-2);
end
for j=5:6
    X(:,j) = X2(:,j-4);
end
for j=7:8
    X(:,j) = X3(:,j-6);
end
for j=9:10
    X(:,j) = X4(:,j-8);
end
for j=11:12
    X(:,j) = X5(:,j-10);
end
for j=13:14
    X(:,j) = X6(:,j-12);
end
for j=15:16
    X(:,j) = X7(:,j-14);
end
for j=17:18
    X(:,j) = X8(:,j-16);
end
for j=19:20
    X(:,j) = X9(:,j-18);
end
for j=21:22
    X(:,j) = X10(:,j-20);
end
for j=23:24
    X(:,j) = X11(:,j-22);
end
for j=25:26
    X(:,j) = X12(:,j-24);
end
for j=27:28
    X(:,j) = X13(:,j-26);
end
for j=29:30
    X(:,j) = X14(:,j-28);
end
for j=31:32
    X(:,j) = X15(:,j-30);
end
for j=33:34
    X(:,j) = X16(:,j-32);
end
for j=35:36
    X(:,j) = X17(:,j-34);
end
for j=37:38
    X(:,j) = X18(:,j-36);
end
for j=39:40
    X(:,j) = X19(:,j-38);
end
for j=41:42
    X(:,j) = X20(:,j-40);
end
for j=43:44
    X(:,j) = X21(:,j-42);
end
for j=45:46
    X(:,j) = X22(:,j-44);
end
for j=47:48
    X(:,j) = X23(:,j-46);
end




XX = zeros(p,n_sample);
simu_type = zeros(p,n_sample);
for i = 1:p
    for j = 1:n_sample
        XX(i,j) = X(i,2*j);
        simu_type(i,j) = X(i,2*j-1);
    end
end
X = XX;

 X1 = zeros(p,2 * n_sample);

for f = 1:n_sample
    for i = 1:p
        X1(i,2*f-1) = simu_type(i,f);
        X1(i,2*f) = X(i,f);
    end
end

%compute mean
X_= X1;
z1 = zeros(p,5);
z2 = zeros(p,5);
z3 = zeros(p,5);
z4 = zeros(p,5);
for i = 1:p
    for j = 1:2:n_training
        if X_(i,j) == 1
            z1(i,1) = z1(i,1) + 1;
            z1(i,2) = z1(i,2) + X_(i,j+1);
        else
            if X_(i,j) == 2
                z2(i,1) = z2(i,1) + 1;
                z2(i,2) = z2(i,2) + X_(i,j+1);
            else
                if X_(i,j) == 3
                    z3(i,1) = z3(i,1) + 1;
                    z3(i,2) = z3(i,2) + X_(i,j+1);
                else
                    if X_(i,j) == 4
                        z4(i,1) = z4(i,1) + 1;
                        z4(i,2) = z4(i,2) + X_(i,j+1);
                    end
                end
            end
        end
    end
end



z1_tot = 0; z2_tot = 0; z3_tot = 0; z4_tot = 0;
z1_count = 0; z2_count = 0; z3_count = 0; z4_count = 0; 
for i = 1:p
    z1_tot = z1_tot + z1(i,2);
    z1_count = z1_count + z1(i,1);
    z2_tot = z2_tot + z2(i,2);
    z2_count = z2_count + z2(i,1);
    z3_tot = z3_tot + z3(i,2);
    z3_count = z3_count + z3(i,1);
    z4_tot = z4_tot + z4(i,2);
    z4_count = z4_count + z4(i,1);
end

z1_mean = z1_tot/z1_count;
z2_mean = z2_tot/z2_count;
z3_mean = z3_tot/z3_count;
z4_mean = z4_tot/z4_count;

z_mean = (z1_mean + z2_mean + z3_mean + z4_mean)/4;

%compute standard error


std_total = zeros(p,1);
for i = 1:p
    for j = 1:2:n_training
         std_total(i,1) = std_total(i,1) + (X_(i,j+1) - z_mean)^2;
    end
    std_total(i,1) = sqrt(std_total(i,1)/n_training);
end



% std_1 = zeros(p,1);
% std_2 = zeros(p,1);
% std_3 = zeros(p,1);
% 
% for i = 1:54
%     for j = 1:2:n_training
%         if X_(i,j) == 1
%             std_1(i,1) = std_1(i,1) + 1;
%             std_1(i,2) = std_1(i,2) + X_(i,j+1);
%         else
%             if X_(i,j) == 2
%                 std_2(i,1) = std_2(i,1) + 1;
%                 z2(i,2) = z2(i,2) + X_(i,j+1);
%             else
%                 if X_(i,j) == 3
%                     z3(i,1) = z3(i,1) + 1;
%                     z3(i,2) = z3(i,2) + X_(i,j+1);
%                 end
%             end
%         end
%     end
% end


%Input distance matrix decided by HEAs structure.

Stru1 = xlsread('Stru_new.xlsx');


beta = zeros(10,n_insample);
cov_s = zeros(p,1);
count_beta_zero = zeros(10,1);   




%Compute single covariance matrix, centralization 
for v = 1:n_insample
    for i = 1:p
        if X_(i,2*v-1) == 1
            cov_s(i,1) = (X(i,v) - z1(i,3))/std_total(i,1);
        end
        if X_(i,2*v-1) == 2
            cov_s(i,1) = (X(i,v) - z2(i,3))/std_total(i,1);
        end
        if X_(i,2*v-1) == 3
            cov_s(i,1) = (X(i,v) - z3(i,3))/std_total(i,1);
        end
        if X_(i,2*v-1) == 4
            cov_s(i,1) = (X(i,v) - z4(i,3))/std_total(i,1);
        end
    end
    
%     for i = 1:54
%         if X_(i,2*v-1) == 1
%             cov_s(i,1) = X(i,v) - z1(i,3);
%         end
%         if X_(i,2*v-1) == 2
%             cov_s(i,1) = X(i,v) - z2(i,3);
%         end
%         if X_(i,2*v-1) == 3
%             cov_s(i,1) = X(i,v) - z3(i,3);
%         end
%     end

    cov = cov_s * cov_s'; %sample covariance matrix
    S = cov;
    X0 = eye(p);
    Y0 = eye(p);
    for i = 1:p
        Stru1(i,i) = 1;
    end
    Xr = Stru1;
    for i = 1:(p-1)
        for j = (i+1):p
            if Xr(i,j) ~= 0
                Y0(i,j) = 0.01;
            end
        end
    end
    for i = 2:p
        for j = 1:(i-1)
            if Xr(i,j) ~= 0
                Y0(i,j) = 0.01;
            end
        end
    end
    for i = 1:p
        Y0(i,i) = 1;
    end
    Ty = X_(:,2 * v - 1);
    Xr = Stru1;
    Dis = ones(p);

    %Compute inverse covariance matrix by CARGO
    nu = p+1;
        c = 1;
        B_prior = c * eye(p);
   [B,T,obj_inner,obj_outer,X_Y] = CARGO(X0, Xr,Y0, S, p,Ty,nu,B_prior);
    
    %Compute beta
    count_beta = zeros(10,1);
    for i = 1:p
        for j = (i+1):p
            if T(i,j) ~= 0
                if X_(i,2*v-1) == 1 && X_(j,2*v-1) == 1
                    beta(1,v) = beta(1,v) - T(i,j)/T(i,i);
                    count_beta(1,1) = count_beta(1,1) + 1;
                end
                if X_(i,2*v-1) == 2 && X_(j,2*v-1) == 2
                    beta(2,v) = beta(2,v) - T(i,j)/T(i,i);
                    count_beta(2,1) = count_beta(2,1) + 1;
                end
                if X_(i,2*v-1) == 3 && X_(j,2*v-1) == 3
                    beta(3,v) = beta(3,v) - T(i,j)/T(i,i);
                    count_beta(3,1) = count_beta(3,1) + 1;
                end
                if X_(i,2*v-1) == 4 && X_(j,2*v-1) == 4
                    beta(4,v) = beta(4,v) - T(i,j)/T(i,i);
                    count_beta(4,1) = count_beta(4,1) + 1;
                end
                if X_(i,2*v-1) == 1 && X_(j,2*v-1) == 2
                    beta(5,v) = beta(5,v) - T(i,j)/T(i,i);
                    count_beta(5,1) = count_beta(5,1) + 1;
                end
                if X_(i,2*v-1) == 2 && X_(j,2*v-1) == 1
                    beta(5,v) = beta(5,v) - T(i,j)/T(i,i);
                    count_beta(5,1) = count_beta(5,1) + 1;
                end
                if X_(i,2*v-1) == 1 && X_(j,2*v-1) == 3
                    beta(6,v) = beta(6,v) - T(i,j)/T(i,i);
                    count_beta(6,1) = count_beta(6,1) + 1;
                end
                if X_(i,2*v-1) == 3 && X_(j,2*v-1) == 1
                    beta(6,v) = beta(6,v) - T(i,j)/T(i,i);
                    count_beta(6,1) = count_beta(6,1) + 1;
                end
                if X_(i,2*v-1) == 1 && X_(j,2*v-1) == 4
                    beta(7,v) = beta(7,v) - T(i,j)/T(i,i);
                    count_beta(7,1) = count_beta(7,1) + 1;
                end
                if X_(i,2*v-1) == 4 && X_(j,2*v-1) == 1
                    beta(7,v) = beta(7,v) - T(i,j)/T(i,i);
                    count_beta(7,1) = count_beta(7,1) + 1;
                end
                if X_(i,2*v-1) == 2 && X_(j,2*v-1) == 3
                    beta(8,v) = beta(8,v) - T(i,j)/T(i,i);
                    count_beta(8,1) = count_beta(8,1) + 1;
                end
                if X_(i,2*v-1) == 3 && X_(j,2*v-1) == 2
                    beta(8,v) = beta(8,v) - T(i,j)/T(i,i);
                    count_beta(8,1) = count_beta(8,1) + 1;
                end
                if X_(i,2*v-1) == 2 && X_(j,2*v-1) == 4
                    beta(9,v) = beta(9,v) - T(i,j)/T(i,i);
                    count_beta(9,1) = count_beta(9,1) + 1;
                end
                if X_(i,2*v-1) == 4 && X_(j,2*v-1) == 2
                    beta(9,v) = beta(9,v) - T(i,j)/T(i,i);
                    count_beta(9,1) = count_beta(9,1) + 1;
                end
                if X_(i,2*v-1) == 3 && X_(j,2*v-1) == 4
                    beta(10,v) = beta(10,v) - T(i,j)/T(i,i);
                    count_beta(10,1) = count_beta(10,1) + 1;
                end
                if X_(i,2*v-1) == 4 && X_(j,2*v-1) == 3
                    beta(10,v) = beta(10,v) - T(i,j)/T(i,i);
                    count_beta(10,1) = count_beta(10,1) + 1;
                end
            end
        end
    end
    for h = 1:10
        if count_beta(h,1) == 0
            count_beta(h,1) = 1;
            count_beta_zero(h,1) = count_beta_zero(h,1) + 1;
        end
    end
    for j = 1:10
        beta(j,v) = beta(j,v)/ count_beta(j,1);
    end
end

beta_aver = zeros(10,1);
beta_sum = zeros(10,1);
for i = 1:n_insample
    for j = 1:10
        beta_sum(j,1) = beta_sum(j,1) + beta(j,i);
    end
end

for i = 1:10
    if count_beta_zero(i,1) == 0
        beta_aver(i,1) = beta_sum(i,1)/n_insample;
    else
        beta_aver(i,1) = beta_sum(i,1)/(n_insample - count_beta_zero(i,1));
    end
end

%plot Iterative step VS ||X-Y||_F and Iterative step VS F(X)
x = linspace(1,35,35);

figure(1)
plot(x,X_Y)
xlabel('Iterative step','FontSize',15)
ylabel('{||X-Y||_F}','FontSize',15);
box off
%hold on
figure(2)
plot(x,obj_outer)
xlabel('Iterative step','FontSize',15)
ylabel('Objective value: F(X)','FontSize',15)
box off



 