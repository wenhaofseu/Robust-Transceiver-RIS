% Please refer to the following reference:
% K. Xu, S. Gong, M. Cui, G. Zhang and S. Ma, 
% "Statistically Robust Transceiver Design for Multi-RIS Assisted Multi-User MIMO Systems," 
% in IEEE Communications Letters, vol. 26, no. 6, pp. 1428-1432, June 2022.
function [max_mse,B,Phi,R,sum_mse,iter_num,CPUtime] = AO_SumMSE_ElementWiseBCDSolver(Hhat_sd,Hhat_sr,Hhat_rd,...
    Sigma_sd,Sigma_sr,Sigma_rd,Psi_sd,Psi_sr,Psi_rd,...
    Pt,sigma_q,Nt,scalar_M,Nr,D,phi,iter_max,epsilon)
    Sigma_sr_sqrt = Sigma_sr^(1/2);Psi_sr_sqrt = Psi_sr^(1/2);
    Sigma_rd_sqrt = Sigma_rd^(1/2);
    Pi1 = Sigma_sr;
    Pi2 = Sigma_sd;
    Xi2 = Psi_sd;
    Xi3 = Psi_rd;
    xi2 = trace(Xi2)/Nr;
    xi3 = trace(Xi3)/Nr;
    B = [eye(D);zeros(Nt-D,D)];%intial
    B = B*Pt/norm(B)^2;
    last_mse = 10;
    tic;
    for iter = 1:iter_max
        Phi = diag(phi);
        Hhat_eq = Hhat_sd+Hhat_rd*Phi*Hhat_sr;
        temp = Sigma_rd_sqrt*Phi*Hhat_sr;
        Pi3 = temp'*temp+trace(Phi*Psi_sr*Phi'*Sigma_rd)*Sigma_sr;
        temp = Hhat_rd*Phi*Psi_sr_sqrt;
        Xi1 = temp*temp';
        xi1 = trace(Xi1)/Nr;
        Upsilonprime = (sigma_q+trace(B*B'*Pi1)*xi1+trace(B*B'*Pi2)*xi2+trace(B*B'*Pi3)*xi3)*eye(Nr);
        temp = Hhat_eq*B;
        R = (temp*temp'+Upsilonprime )\temp;
        temp = Hhat_eq'*R;
        C = trace(R*R')*(xi1*Pi1+xi2*Pi2+xi3*Pi3)+temp*temp';
        [Uc,LAMBDAc]=eig(C);
        A = Uc'*(temp*temp')*Uc;
        %% check lambda = 0 is optimal?
        flag = 0;
        if rank(C)==Nt
            temp1 = 0;
            for n = 1:Nt
                temp1 = temp1+A(n,n)/LAMBDAc(n,n)^2;
            end
            if temp1<=Pt
                flag = 1;
            end
        end
        %% calculate lambda
        lambda_min = 0;
        lambda_max = sqrt(trace(A)/Pt);
        while(flag==0)
            lambda = (lambda_min+lambda_max)/2;
            temp2 = 0;
            for n = 1:Nt
                temp2 = temp2+A(n,n)/(LAMBDAc(n,n)+lambda)^2;
            end
            temp2 = real(temp2);
            if(temp2<Pt)
                lambda_max = lambda;
            else
                lambda_min = lambda;
            end
            %% bisection converage
            if abs(lambda_max-lambda_min)<1e-3 && abs(temp2-Pt)<1e-5
                break
            elseif abs(lambda_max-lambda_min)<1e-7 && lambda_min==0
                break
            end
        end
        B = (C+lambda*eye(Nt))\temp;
        Upsilonprime = (sigma_q+trace(B*B'*Pi1)*xi1+trace(B*B'*Pi2)*xi2+trace(B*B'*Pi3)*xi3)*eye(Nr);
        temp = R'*Hhat_eq*B;
        Tprime = eye(D)+temp*temp'-temp-temp'+R'*Upsilonprime*R;
        mse = real(trace(Tprime));
        if abs(mse-last_mse)/mse<=epsilon
            sum_mse = mse;
            iter_num = iter;
            break
        end
        last_mse = mse;
        Gamma1 = KR(R'*Hhat_rd,(Hhat_sr*B).');
        Gamma2 = norm(R,'fro')*norm(Sigma_sr_sqrt*B,'fro')/sqrt(Nr)...
            *KR(Hhat_rd,Psi_sr_sqrt.');
        Gamma3 = sqrt(xi3)*norm(R,'fro')*KR(Sigma_rd_sqrt,(Hhat_sr*B).');
        Gamma4 = sqrt(xi3)*norm(R,'fro')*norm(Sigma_sr_sqrt*B,'fro')...
            *KR(Sigma_rd_sqrt,Psi_sr_sqrt.');
        Gamma = Gamma1'*Gamma1+Gamma2'*Gamma2...
            +Gamma3'*Gamma3+Gamma4'*Gamma4;
        gamma = Gamma1'*vec((R'*Hhat_sd*B-eye(D)).');
        phi = ElementWiseBCDupdate(Gamma,gamma,phi,scalar_M);
    end
    CPUtime = toc;
    Upsilon = sigma_q*eye(Nr)...
        +trace(B*B'*Pi1)*Xi1+trace(B*B'*Pi2)*Xi2+trace(B*B'*Pi3)*Xi3;
    temp = Hhat_eq*B;
    R = (temp*temp'+Upsilon)\temp;
    T = eye(D)-temp'*R;
    [V2,~] = eig(T);
    B = B*V2;
    R = R*V2;
    max_mse = max(real(diag(T)));
end

