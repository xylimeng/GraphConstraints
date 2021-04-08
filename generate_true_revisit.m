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
% Ty:      true atom type vector
function[Xr,Xrt,Ty] = generate_true(p,sparsity,a)
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
Sq = zeros(p,p);
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
                C(i,j) =  -0.2;
            end
            if Ty(i) == 2 && Ty(j) == 2
                C(i,j) =  -0.4;
            end
            if Ty(i) == 3 && Ty(j) == 3
                C(i,j) =  0.4;
            end
            if Ty(i) == 1 && Ty(j) == 2
                C(i,j) =  -0.42;
            end
            if Ty(i) == 2 && Ty(j) == 1
                C(i,j) =  -0.42;
            end
            if Ty(i) == 1 && Ty(j) == 3
                C(i,j) =  0.25;
            end
            if Ty(i) == 3 && Ty(j) == 1
                C(i,j) =  0.25;
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
%                 Sq(i,j) =  0.01;
%             end
%             if Ty(i) == 2 && Ty(j) == 2
%                 Sq(i,j) =  0.01;
%             end
%             if Ty(i) == 3 && Ty(j) == 3
%                 Sq(i,j) =  0.01;
%             end
%             if Ty(i) == 1 && Ty(j) == 2
%                 Sq(i,j) =  0.01;
%             end
%             if Ty(i) == 2 && Ty(j) == 1
%                 Sq(i,j) =  0.01;
%             end
%             if Ty(i) == 1 && Ty(j) == 3
%                 Sq(i,j) =  0.01;
%             end
%             if Ty(i) == 3 && Ty(j) == 1
%                 Sq(i,j) =  0.01;
%             end
%             if Ty(i) == 2 && Ty(j) == 3
%                 Sq(i,j) =  0.01;
%             end
%             if Ty(i) == 3 && Ty(j) == 2
%                 Sq(i,j) =  0.01;
%             end
%         end
%     end
% end
% for i = 1:p
%     Sq(i,i) = 1;%4;
% end


for i = 2:p
    for j = 1:i-1
        C(i,j) = C(j,i);
    end
end

for i = 2:p
    for j = 1:i-1
        Sq(i,j) = Sq(j,i);
    end
end
while min(eig(C)) <= 0
    for i = 1:p
        C(i,i) = C(i,i) + 0.1;
    end
end
Xr = C;
Xrt = inv(Xr);
 
