% Please refer to the following literature:
% D. P. Palomar, J. M. Cioffi and M. A. Lagunas, 
% "Joint Tx-Rx beamforming design for multicarrier MIMO channels: a unified framework for convex optimization," 
% in IEEE Transactions on Signal Processing, vol. 51, no. 9, pp. 2381-2401, Sept. 2003.
function [p] = waterfilling(Pt,eigen,D)
    %Pt is total transmit power
    %eigen is Nt-dimensional vector 
    %L = min(D,rank(H))
    %p is allocated power and d-dimensional vector
    r = length(find(eigen>0));
    %delet zero element
    eigen(1:length(eigen)-r) = [];
    re_eigen = [zeros(max(D-r,0),1);eigen(1+max(r-D,0):r)];
    re_eigen1 = 1./re_eigen(D+1-min(D,r):D);
    p(1:max(D-r,0))=0;
    p(1+max(D-r,0):D) = min_max_wf(Pt,re_eigen1);    
end

