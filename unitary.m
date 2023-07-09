function [U] = unitary(lamda,label,order)
    %   order = matrix order
    %   label is a vector,which has two elements
    I = eye(order);
    Q = perm(I,label(1),label(2));
    T = lamda*I+(1-lamda)*Q;
    U = zeros(order);
    for i=1:order
        for j=1:order
            if i<=j
                U(i,j)=sqrt(T(i,j));
            else
                U(i,j)=-1*sqrt(T(i,j));
            end
        end
    end
end

