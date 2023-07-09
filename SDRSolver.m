function [max_mse,B,Phi,R] = SDRSolver(Hhat_sd,Hhat_sr,Hhat_rd,...
    Sigma_sd,Sigma_sr,Sigma_rd,Psi_sd,Psi_sr,Psi_rd,...
    Pt,sigma_q,Nt,M,Nr,D,phi,iter_max,epsilon,GaussianRandomizationNumber)
    Sigma_sd_sqrt = Sigma_sd^(1/2);Psi_sd_sqrt = Psi_sd^(1/2);
    Sigma_sr_sqrt = Sigma_sr^(1/2);Psi_sr_sqrt = Psi_sr^(1/2);
    Sigma_rd_sqrt = Sigma_rd^(1/2);Psi_rd_sqrt = Psi_rd^(1/2);
    Pi1 = Sigma_sr;
    Pi2 = Sigma_sd;
    Xi2 = Psi_sd;
    Xi3 = Psi_rd;
    mat2 = kron(eye(D),Sigma_sr_sqrt);
    mat3 = kron(eye(D),Sigma_sd_sqrt);
    I_D = eye(D);
    Phi = diag(phi);
	%% Use phi to initialize B
    Hhat_eq = Hhat_sd+Hhat_rd*Phi*Hhat_sr;
    [U,Sigma] = eig(Hhat_eq'*Hhat_eq/sigma_q);
    rank_Hhat_eq = rank(Hhat_eq'*Hhat_eq/sigma_q);
    L = min(D,rank_Hhat_eq);
    U1 = U(:,Nt-min(D,L)+1:Nt);
    g = waterfilling(Pt,diag(Sigma),D);
    GG = diag(sqrt(g));%   noting sqrt
    Sigma1 = [zeros(L,D-L),GG(1+max(D-L,0):max(D,L),1+max(D-L,0):max(D,L))];
    B = U1*Sigma1;
    T = inv(eye(D)+(Hhat_eq*B)'*(Hhat_eq*B)/sigma_q);
    V = rotation(diag(T));
    B = B*V;
    last_mse = 10; %record sum-mse
    for iter = 1:iter_max
		%% Update R
        Hhat_eq = Hhat_sd+Hhat_rd*Phi*Hhat_sr;
        temp = Sigma_rd_sqrt*Phi*Hhat_sr;
        Pi3 = temp'*temp+trace(Phi*Psi_sr*Phi'*Sigma_rd)*Sigma_sr;
        temp = Hhat_rd*Phi*Psi_sr_sqrt;
        Xi1 = temp*temp';
        Upsilon = sigma_q*eye(Nr)...
        +trace(B*B'*Pi1)*Xi1+trace(B*B'*Pi2)*Xi2+trace(B*B'*Pi3)*Xi3;
        temp = Hhat_eq*B;
        R = (temp*temp'+Upsilon)\temp;
        R_H = R';
		%% Update B (SOCP)
        mat4 = kron(eye(D),Sigma_rd_sqrt*Phi*Hhat_sr);
        cvx_begin quiet
        variable B(Nt, D) complex;
        variable tau1;
        expression m(D+2*Nt*D+M*D+Nr,D);
        for d = 1:D
            r_d = R(:,d);
            %noting:mat1,K1,K2,K3,K4 is related to R_H(i,:)
            mat1 = kron(I_D,R_H(d,:)*Hhat_eq); 
            k1=norm(Psi_sr_sqrt*Phi'*Hhat_rd'*r_d,2)^2;
            k3=norm(Psi_sd_sqrt*r_d,2)^2;
            k4=norm(Psi_rd_sqrt*r_d,2)^2;
            k2=trace(Phi*Psi_sr*Phi'*Sigma_rd)*k4;
            %describe SOC constraints
            m(:,d)=[...
                mat1*vec(B)-I_D(:,d);...
                sqrt(k1+k2)*mat2*vec(B);...
                sqrt(k3)*mat3*vec(B);...
                sqrt(k4)*mat4*vec(B);...
                sqrt(sigma_q)*r_d...
                ];
        end
        minimize(tau1);
        subject to
        for d = 1:D
            norm(m(:,d))<=tau1;%noting sqrt
        end
        norm(B,'fro')<=sqrt(Pt);
        cvx_end
        mse = tau1^2;
		%% Convergence criteria
        if abs(mse-last_mse)/mse<=epsilon
            break
        end
		%% store variable max-mse
        last_mse = mse;
        
        %% Update Phi (SDR)
		Z=zeros(M+2,M+2);
        for d = 1:D
            r_d = R(:,d);
            Gamma1 = KR(R_H(d,:)*Hhat_rd,(Hhat_sr*B).');
            Gamma2 = norm(Sigma_sr_sqrt*B,'fro')...
                *KR(R_H(d,:)*Hhat_rd,Psi_sr_sqrt.');
            Gamma3 = norm(Psi_rd_sqrt*r_d)*KR(Sigma_rd_sqrt,(Hhat_sr*B).');
            Gamma4 = norm(Psi_rd_sqrt*r_d)*norm(Sigma_sr_sqrt*B,'fro')...
                *KR(Sigma_rd_sqrt,Psi_sr_sqrt.');
            Gamma = Gamma1'*Gamma1+Gamma2'*Gamma2+Gamma3'*Gamma3+Gamma4'*Gamma4;
            temp = R_H(d,:)*Hhat_sd*B-I_D(:,d)';
            gamma = Gamma1'*(temp.');
            scalar_gamma = norm(temp)^2+sigma_q*norm(r_d)^2+...
                norm(Psi_sd_sqrt*r_d)^2*norm(Sigma_sd_sqrt*B,'fro')^2;
            Z(:,:,d)=[...
                Gamma,       zeros(M,1),        gamma;...
                zeros(1,M),     0,               -0.5;...
                gamma',        -0.5,      scalar_gamma...
                ];
        end
        cvx_begin quiet
        variable X(M+2,M+2) complex semidefinite;
        variable tau;
        minimize(tau);
        subject to
        for d = 1:D
            real(trace(X*Z(:,:,d))) <= 0;
        end
        diag(X(1:M,1:M))==1;
        X(M+2,M+2)==1;
        X(M+1,M+1)==tau;
        tau>=0;
        cvx_end
		
        %%  Gaussian randomization
        [U,Sigma]=eig(X);
        phi_opt=zeros(M,GaussianRandomizationNumber);mse=zeros(D,GaussianRandomizationNumber);
        mse_max=zeros(GaussianRandomizationNumber,1);
        for num=1:GaussianRandomizationNumber
            w1 = sqrt(1/2)*(randn(M+2,1)+1j*randn(M+2,1));
            w2 = U*sqrt(Sigma)*w1;
            phi_opt(:,num) = exp(1j*angle(w2(1:M)/w2(M+2)));
            temp = [phi_opt(:,num);0;1];
            temp = temp*temp';
            for d = 1:D
                mse(d,num)=trace(temp*Z(:,:,d));
            end
            mse_max(num)=max(mse(:,num),[],1); 
        end
		%Select the one with the smallest MAX-MSE.
        [~,idx]=min(mse_max);
        phi = phi_opt(:,idx);   %opt phi
        Phi = diag(phi);
    end
	%% record per-stream MAX-MSE
    max_mse = mse;
end

