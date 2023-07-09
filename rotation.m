% s.t. U'*A*U=B ,A is diag matrix and diag(B) is equal
% Please refer to the following literature:
% P. Viswanath and V. Anantharam, "Optimal sequences and sum capacity of synchronous CDMA systems," 
% in IEEE Transactions on Information Theory, vol. 45, no. 6, pp. 1984-1991, Sept. 1999.
function [U] = rotation(a)
    order = length(a);
    U = eye(order);
    aver = sum(a)/order;
    above_index = find(a>aver);
    below_index = find(a<aver);
    equal_index = find(a==aver);
    for i=1:order-1
        %calculate lamda(i)
        %s.t. lamda(i)*a(above_index(1))+(1-lamda(i))*a(below_index(1))=aver
        lamda(i)=(aver-a(below_index(1)))/(a(above_index(1))-a(below_index(1)));
        label(:,i) = [above_index(1);below_index(1)];
        a(below_index(1))=a(above_index(1))+a(below_index(1))-aver;
        a(above_index(1))=aver;
        %%  update
        above_index = find(a>aver);
        below_index = find(a<aver);
        equal_index = find(a==aver);
        if length(equal_index)==order
            break
        end
    end
    num = i;
    for j = 1:num
        U = U*unitary(lamda(j),label(:,j),order);
    end
end

