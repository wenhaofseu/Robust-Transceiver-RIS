function [max_mse,B,Phi,R] = ManoptSolver(Hhat_sd,Hhat_sr,Hhat_rd,...
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
    last_mse = 10; %record sum-mse
    for iter = 1:iter_max
        Phi = diag(phi);
        temp = Sigma_rd_sqrt*Phi*Hhat_sr;
        Pi3 = temp'*temp+trace(Phi*Psi_sr*Phi'*Sigma_rd)*Sigma_sr;
        temp = Hhat_rd*Phi*Psi_sr_sqrt;
        Xi1 = temp*temp';
        xi1 = trace(Xi1)/Nr;
		%% derive a closed-form solution for the joint B and R design
        M = sigma_q*eye(Nt)+Pt*(xi1*Pi1+xi2*Pi2+xi3*Pi3);
        Inv_M = M^(-1);
        Inv_M_sqrt = M^(-1/2);
        Hhat_eq = Hhat_sd+Hhat_rd*Phi*Hhat_sr;
        temp = Hhat_eq*Inv_M_sqrt;
        [U,Lambda,V] = svd(temp');
        U1 = U(:,1:D);
        V1 = V(:,1:D);
        Lambda1 = Lambda(1:D,1:D);
        for L = D:-1:1
            Lambda1_L = Lambda1(1:L,1:L);
            lambda1_L = real(Lambda1(L,L));
            Inv_Lambda1_L = Lambda1_L^(-1);
            U1_L = U1(:,1:L);
            a1 = trace(Inv_Lambda1_L^2);
            a2 = trace(Inv_Lambda1_L);
            a3 = trace(Inv_Lambda1_L*U1_L'*Inv_M*U1_L);
            a4 = trace(Inv_Lambda1_L^2*U1_L'*Inv_M*U1_L);
            eta = (a2*Pt)/(a3*(Pt+a1)-a2*a4);
            eta = real(eta);
            lambda = (sigma_q*a2*(a3*(Pt+a1)-a2*a4))/(Pt*(Pt+a1)^2);
            lambda = real(lambda);
            if (lambda<=lambda1_L^2*sigma_q/eta)
                break;
            end
        end
        Lambda_B_L = (sqrt(eta*sigma_q/lambda)*Inv_Lambda1_L-eta*Inv_Lambda1_L^2)^(1/2);
        Lambda_R_L = sqrt(lambda/(eta*sigma_q))*Lambda_B_L*Inv_Lambda1_L;
        Lambda_B = [Lambda_B_L,zeros(L,D-L);zeros(D-L,L),zeros(D-L,D-L)];
        Lambda_R = [Lambda_R_L,zeros(L,D-L);zeros(D-L,L),zeros(D-L,D-L)];
        B = Inv_M_sqrt*U1*Lambda_B;
        R = V1*Lambda1*Lambda_R;
%         Upsilon_prime = eta*eye(Nr);
%         temp = R'*Hhat_eq*B;
%         Tprime = eye(D)+temp*temp'-temp-temp'+R'*Upsilon_prime*R;
		%% record sum-mse
        Tprime = inv(eye(D)+Lambda_B'*Lambda1^2*Lambda_B/eta);
        mse = real(trace(Tprime));
		%% Convergence criteria
        if abs(mse-last_mse)/mse<=epsilon
            break
        end
		%% store variable sum-mse
        last_mse = mse;
		%% manifold-based conjugate gradient method
        Gamma1 = KR(R'*Hhat_rd,(Hhat_sr*B).');
        Gamma2 = norm(R,'fro')*norm(Sigma_sr_sqrt*B,'fro')/sqrt(Nr)...
            *KR(Hhat_rd,Psi_sr_sqrt.');
        Gamma3 = sqrt(xi3)*norm(R,'fro')*KR(Sigma_rd_sqrt,(Hhat_sr*B).');
        Gamma4 = sqrt(xi3)*norm(R,'fro')*norm(Sigma_sr_sqrt*B,'fro')...
            *KR(Sigma_rd_sqrt,Psi_sr_sqrt.');
        Gamma = Gamma1'*Gamma1+Gamma2'*Gamma2...
            +Gamma3'*Gamma3+Gamma4'*Gamma4;
        gamma = Gamma1'*vec((R'*Hhat_sd*B-eye(D)).');
        func = @(phi) real(phi'*Gamma*phi)+2*real(gamma'*phi); %objective function related to variable phi
        manopt = complexcirclefactory(scalar_M);
        problem.M = manopt;
        problem.cost = func;
        problem.egrad = @(phi) 2*Gamma*phi+2*gamma; %Multiply the partial derivative by the coefficient 2 to obtain the gradient
        options.verbosity = 0;
        [phi, ~, ~, ~] = conjugategradient(problem,phi,options);
    end
	%% Use the true values of matrix Upsilon instead of its approximations.
    Upsilon = sigma_q*eye(Nr)...
        +trace(B*B'*Pi1)*Xi1+trace(B*B'*Pi2)*Xi2+trace(B*B'*Pi3)*Xi3;
	%% fix B, then update R and record MSE matrix
    temp = Hhat_eq*B; 
    R = (temp*temp'+Upsilon)\temp;
    T = eye(D)-temp'*R;
	%% Diagonalize the MSE matrix
    [V2,~] = eig(T);
    B = B*V2;
    R = R*V2;
    Upsilon = sigma_q*eye(Nr)...
        +trace(B*B'*Pi1)*Xi1+trace(B*B'*Pi2)*Xi2+trace(B*B'*Pi3)*Xi3;
    temp = R'*Hhat_eq*B;
    T = eye(D)+temp*temp'-temp-temp'+R'*Upsilon*R;
	%% Make the main diagonal elements of the MSE matrix equal
    V2 = rotation(diag(T)); %
    B = B*V2;
    R = R*V2;
	%% record per-stream MAX-MSE
    max_mse = mean(real(diag(T)));
end

