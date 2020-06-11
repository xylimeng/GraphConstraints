% function [S,X0,Y0,c1,h1] = generate_samples(h,p,Xrt,samplenumber,a1,b1)
%
% input:
%
%a:select a seed for random generating samples based on Xr
%samplenumber: number of samples
%
% output:
%
% S: the correlation matrix
% X0:initial set-up of variable X;
% Y0:initial set-up of variable Y;
% r_samples: the generated samples;
%
function [r_samples,S,X0,Y0] = generatesamples(h,p,Xr,Xrt,samplenumber)
randn('seed',h);
m = zeros(p,1)';
r = mvnrnd(m,Xrt,samplenumber); 
r_samples = r';
%csvwrite('simu_samples.csv',r_samples);
%S = cov(r);      
%S = corrcoef(r');   
S = r'*r;   %单个样本下的协方差计算
%S = randn(p)
X0 = eye(p);
Y0 = eye(p);
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
     Y0(i,i) = Xr(i,i);
end

 
% if mod(h,a1) == 0     
%     h1 = fix(h/a1);
% else
%     h1 = fix(h/a1) + 1;
% end