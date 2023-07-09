%%   fuction:generate Permutation matrix 
function [A] = perm(A,m,n)
    temp = A(:,m);
    A(:,m) = A(:,n);
    A(:,n) = temp;
end

