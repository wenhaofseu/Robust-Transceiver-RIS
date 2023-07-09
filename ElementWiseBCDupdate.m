% Please refer to the following reference:
% K. Xu, S. Gong, M. Cui, G. Zhang and S. Ma, 
% "Statistically Robust Transceiver Design for Multi-RIS Assisted Multi-User MIMO Systems," 
% in IEEE Communications Letters, vol. 26, no. 6, pp. 1428-1432, June 2022.(See Section III B)
%min x'*A*x+b'*x+x'*b
%s.t. \lvert xm \rvert = 1, m=1,2,\cdots,M
%fixed {xi}_{i=1,i\neq m}, optimitize xm
%min 2\Re{conj(xm)*(\sum_{i=1,i\neq m}Ami*xi+bm)}
function x = ElementWiseBCDupdate(A,b,x,M)
    for m = 1:M
        %update x(m)
        x(m) = -exp(1i*angle(A(m,:)*x-A(m,m)*x(m)+b(m)));
    end
end

