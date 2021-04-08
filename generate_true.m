% generate true inverse covariance matrix
% function[Xr,Xrt,S1] = generate_true(p,sparsity,a)
%
% input:
%
% p: dimension of true inverse covariance matrix
% sparsity: fraction of nonzero elements in the inverse covariance matrix
% a: select a seed for random generating a true inverse covariance matrix
%
% output:
%
% Xr :     true inverse covariance matrix
% Xrt:     true covariance matrix
% S1:      number of non-zero elements of Xrt
% Ty:      true atom type vector
function[Xr,Xrt,S1,Ty] = generate_true(p,sparsity,a)
rand('seed',a);
Ty = zeros(p,1);
for i = 1:p
    Ty(i) = rand(1);
    if Ty(i) <= 1/3
        Ty(i) = 1;
    end
    if 1/3 < Ty(i) && Ty(i) <= 2/3
        Ty(i) = 2;
    end
    if 2/3 < Ty(i) && Ty(i) < 1
        Ty(i) = 3;
    end
end

A = sprand(p,p,sparsity);
U = zeros(p,p);
C = zeros(p,p);
for i = 1:p
    for j = 1:p
        U(i,j) = A(i,j);
    end
end
for i = 1:p
    for j = 1:p
        if U(i,j) >= 0.5
            U(i,j) = 1;
        end
    end
end
for i = 1:p
    for j = 1:p
        if U(i,j) ~= 0 && U(i,j) <= 0.5
            U(i,j) = -1;
        end
    end
end

%U = xlsread('true_structure.xlsx');

for i = 1:p
    for j = i:p
        if U(i,j) ~= 0
            if Ty(i) == 1 && Ty(j) == 1
                C(i,j) =  0.2;
            end
            if Ty(i) == 2 && Ty(j) == 2
                C(i,j) =  0.4;
            end
            if Ty(i) == 3 && Ty(j) == 3
                C(i,j) =  0.8;
            end
            if Ty(i) == 1 && Ty(j) == 2
                C(i,j) =  -0.42;
            end
            if Ty(i) == 2 && Ty(j) == 1
                C(i,j) =  -0.42;
            end
            if Ty(i) == 1 && Ty(j) == 3
                C(i,j) =  -0.65;
            end
            if Ty(i) == 3 && Ty(j) == 1
                C(i,j) =  -0.65;
            end
            if Ty(i) == 2 && Ty(j) == 3
                C(i,j) =  0.5;
            end
            if Ty(i) == 3 && Ty(j) == 2
                C(i,j) =  0.5;
            end
        end
    end
end
for i = 1:p
    C(i,i) = 4;%4;
end



% for i = 1:p
%     for j = i:p
%         if U(i,j) ~= 0
%             if Ty(i) == 1 && Ty(j) == 1
%                 C(i,j) =  0.1;
%             end
%             if Ty(i) == 2 && Ty(j) == 2
%                 C(i,j) =  0.1;
%             end
%             if Ty(i) == 3 && Ty(j) == 3
%                 C(i,j) =  0.1;
%             end
%             if Ty(i) == 1 && Ty(j) == 2
%                 C(i,j) =  -0.08;
%             end
%             if Ty(i) == 2 && Ty(j) == 1
%                 C(i,j) =  -0.08;
%             end
%             if Ty(i) == 1 && Ty(j) == 3
%                 C(i,j) =  -0.1;
%             end
%             if Ty(i) == 3 && Ty(j) == 1
%                 C(i,j) =  -0.1;
%             end
%             if Ty(i) == 2 && Ty(j) == 3
%                 C(i,j) =  0.2;
%             end
%             if Ty(i) == 3 && Ty(j) == 2
%                 C(i,j) =  0.2;
%             end
%         end
%     end
% end
% for i = 1:p
%     C(i,i) = 1;%4;
% end


for i = 2:p
    for j = 1:i-1
        C(i,j) = C(j,i);
    end
end
while min(eig(C)) <= 0
    for i = 1:p
        C(i,i) = C(i,i) + 0.1;
    end
end
Xr = C;
%Xr = xlsread('true_structure.xlsx');
Xrt = inv(Xr);
%csvwrite('Type.csv',Ty);
%csvwrite('Xr.csv',Xr);       %save Xr in "Xr.csv"  
%csvwrite('Xrt.csv',Xrt);        %save Xrt in "Xrt.csv"
S1 = sum(sum(Xr~=0,2));      %sparsity rate

