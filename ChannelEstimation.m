clc;clear;close all;
rng(0);
Par = load_parameters();
load('RISchannel.mat','Hsd_all','Hsr_all','Hrd_all',...
    'Hbar_sd','Hbar_sr','Hbar_rd',...
    'R_Tx_sd','R_Rx_sd','R_Tx_sr','R_Rx_sr','R_Tx_rd','R_Rx_rd',...
    'pl_sd','pl_sr','pl_rd');
R_Tx_sd_sqrt = R_Tx_sd^(1/2); R_Rx_sd_sqrt = R_Rx_sd^(1/2);
R_Tx_sr_sqrt = R_Tx_sr^(1/2); R_Rx_sr_sqrt = R_Rx_sr^(1/2);
R_Tx_rd_sqrt = R_Tx_rd^(1/2); R_Rx_rd_sqrt = R_Rx_rd^(1/2);
Sigma_sd = R_Tx_sd;
Sigma_sr = R_Tx_sr;
Sigma_rd = R_Tx_rd;
%% setting number of antennas
Nt = Par.Nt; 
M = Par.M;
Nr = Par.Nr;
%% calculate Inv
Inv_R_Rx_sd = inv(R_Rx_sd);
Inv_R_Rx_sr = inv(R_Rx_sr);
Inv_R_Rx_rd = inv(R_Rx_rd);
%% setting Rician factor
kappa_sd = Par.kappa_sd;
kappa_sr = Par.kappa_sr;
kappa_rd = Par.kappa_rd;
%% calculate radio of LOS and NLOS
kappa2_sd = 1/(1+kappa_sd); kappa1_sd = 1-kappa2_sd;
kappa1_sd = sqrt(kappa1_sd); kappa2_sd = sqrt(kappa2_sd);
kappa2_sr = 1/(1+kappa_sr); kappa1_sr = 1-kappa2_sr;
kappa1_sr = sqrt(kappa1_sr); kappa2_sr = sqrt(kappa2_sr);
kappa2_rd = 1/(1+kappa_rd); kappa1_rd = 1-kappa2_rd;
kappa1_rd = sqrt(kappa1_rd); kappa2_rd = sqrt(kappa2_rd);
%% noise power
sigma_q = Par.sigma_q;
%% pilot power for Hsd
SNR_range = Par.SNR_range;
%% normalized pilot matrix
U0_sd = dftmtx(Nt)/sqrt(Nt);
X0_sd = R_Tx_sd^(-1/2)*U0_sd;
temp_sd = norm(X0_sd,'fro');
X0_sd = X0_sd/temp_sd;
Inv_Usd = temp_sd*U0_sd';
U0_sr = dftmtx(Nt)/sqrt(Nt);
X0_sr = R_Tx_sr^(-1/2)*U0_sr;
temp_sr = norm(X0_sr,'fro');
Inv_Usr = temp_sd*U0_sr';
X0_sr = X0_sr/temp_sr;
U0_rd = dftmtx(M)/sqrt(M);
X0_rd = R_Tx_rd^(-1/2)*U0_rd;
temp_rd = norm(X0_rd,'fro');
X0_rd = X0_rd/temp_rd;
Inv_Urd = temp_rd*U0_rd';
%% Monte Carlo
number = Par.ChannelEstimationNumber;
%% store result
err_Hsd = zeros(length(SNR_range),number);
err_Hsr = zeros(length(SNR_range),number);
err_Hrd = zeros(length(SNR_range),number);
MSE_Hsd_theory = zeros(length(SNR_range),1);
MSE_Hsr_theory = zeros(length(SNR_range),1);
MSE_Hrd_theory = zeros(length(SNR_range),1);
for i = 1:length(SNR_range)
    %% pilot power
    SNR = SNR_range(i);
    P_sd = sigma_q/pl_sd*10.^(SNR/10);
    P_sr = sigma_q/pl_sr*10.^(SNR/10);
    P_rd = sigma_q/pl_rd*10.^(SNR/10);
    sigma_ce_q_sd = temp_sd^2*sigma_q/P_sd;
    sigma_ce_q_sr = temp_sr^2*sigma_q/P_sr;
    sigma_ce_q_rd = temp_rd^2*sigma_q/P_rd;
    sigma_0_q_sd = sigma_ce_q_sd*(1+kappa_sd)/pl_sd;
    sigma_0_q_sr = sigma_ce_q_sr*(1+kappa_sr)/pl_sr;
    sigma_0_q_rd = sigma_ce_q_rd*(1+kappa_rd)/pl_rd;
    Psi_sd = sigma_ce_q_sd.*inv(eye(Nr)+sigma_0_q_sd.*Inv_R_Rx_sd);
    Psi_sr = sigma_ce_q_sr.*inv(eye(M)+sigma_0_q_sr.*Inv_R_Rx_sr);
    Psi_rd = sigma_ce_q_rd.*inv(eye(Nr)+sigma_0_q_rd.*Inv_R_Rx_rd);
    MSE_Hsd_theory(i) = trace(Psi_sd)*trace(Sigma_sd);
    MSE_Hsr_theory(i) = trace(Psi_sr)*trace(Sigma_sr);
    MSE_Hrd_theory(i) = trace(Psi_rd)*trace(Sigma_rd);
    %% transmitted pilot matrix
    Xsd = sqrt(P_sd)*X0_sd; 
    Xsr = sqrt(P_sr)*X0_sr; 
    Xrd = sqrt(P_rd)*X0_rd;
    for num = 1:number
        num
        Hsd = Hsd_all(:,:,num);
        Hsr = Hsr_all(:,:,num);
        Hrd = Hrd_all(:,:,num);
        Nsd = sqrt(sigma_q/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
        Nsr = sqrt(sigma_q/2)*(randn(M,Nt)+1i*randn(M,Nt));
        Nrd = sqrt(sigma_q/2)*(randn(Nr,M)+1i*randn(Nr,M));
        %% received pilot signal
        Ysd = Hsd*Xsd+Nsd;
        Ysr = Hsr*Xsr+Nsr;
        Yrd = Hrd*Xrd+Nrd;
        %% linear observation model
        Htilde_sd_ob = sqrt((1+kappa_sd)/pl_sd)*R_Rx_sd^(-1/2)*Ysd*Inv_Usd/sqrt(P_sd)...
            -sqrt(kappa_sd)*R_Rx_sd^(-1/2)*Hbar_sd*Xsd*Inv_Usd/sqrt(P_sd);
        Htilde_sr_ob = sqrt((1+kappa_sr)/pl_sr)*R_Rx_sr^(-1/2)*Ysr*Inv_Usr/sqrt(P_sr)...
            -sqrt(kappa_sr)*R_Rx_sr^(-1/2)*Hbar_sr*Xsr*Inv_Usr/sqrt(P_sr);
        Htilde_rd_ob = sqrt((1+kappa_rd)/pl_rd)*R_Rx_rd^(-1/2)*Yrd*Inv_Urd/sqrt(P_rd)...
            -sqrt(kappa_rd)*R_Rx_rd^(-1/2)*Hbar_rd*Xrd*Inv_Urd/sqrt(P_rd);
        %% LMMSE estimator of Htilde
        Htilde_sd_hat = (eye(Nr)+sigma_0_q_sd.*Inv_R_Rx_sd)\Htilde_sd_ob;
        Htilde_sr_hat = (eye(M)+sigma_0_q_sr.*Inv_R_Rx_sr)\Htilde_sr_ob;
        Htilde_rd_hat = (eye(Nr)+sigma_0_q_rd.*Inv_R_Rx_rd)\Htilde_rd_ob;
        %% LMMSE estimator of H
        Hhat_sd = sqrt(pl_sd)*(kappa1_sd*Hbar_sd+kappa2_sd*R_Rx_sd_sqrt*Htilde_sd_hat*R_Tx_sd_sqrt);
        Hhat_sr = sqrt(pl_sr)*(kappa1_sr*Hbar_sr+kappa2_sr*R_Rx_sr_sqrt*Htilde_sr_hat*R_Tx_sr_sqrt);
        Hhat_rd = sqrt(pl_rd)*(kappa1_rd*Hbar_rd+kappa2_rd*R_Rx_rd_sqrt*Htilde_rd_hat*R_Tx_rd_sqrt);
        %% calculate error
        err_Hsd(i,num) = norm(Hsd-Hhat_sd,'fro')^2;
        err_Hsr(i,num) = norm(Hsr-Hhat_sr,'fro')^2;
        err_Hrd(i,num) = norm(Hrd-Hhat_rd,'fro')^2;
    end
end
NMSE_Hsd = mean(err_Hsd,2)/(Nr*Nt*pl_sd);
NMSE_Hsr = mean(err_Hsr,2)/(M*Nt*pl_sr);
NMSE_Hrd = mean(err_Hrd,2)/(Nr*M*pl_rd);
NMSE_Hsd_theory = MSE_Hsd_theory/(Nr*Nt*pl_sd);
NMSE_Hsr_theory = MSE_Hsr_theory/(M*Nt*pl_sr);
NMSE_Hrd_theory = MSE_Hrd_theory/(Nr*M*pl_rd);
plot(SNR_range,10*log10(NMSE_Hsd),'rd','linewidth', 2,'Markersize',8);hold on;
plot(SNR_range,10*log10(NMSE_Hsd_theory),'r-','linewidth', 2,'Markersize',8);hold on;
plot(SNR_range,10*log10(NMSE_Hsr),'g<','linewidth', 2,'Markersize',8);hold on;
plot(SNR_range,10*log10(NMSE_Hsr_theory),'g:','linewidth', 2,'Markersize',8);hold on;
plot(SNR_range,10*log10(NMSE_Hrd),'bs','linewidth', 2,'Markersize',8);hold on;
plot(SNR_range,10*log10(NMSE_Hrd_theory),'b--','linewidth', 2,'Markersize',8);hold on;
legend('BS-UE channel, sim.','BS-UE channel, theory',...
    'BS-RIS channel, sim.','BS-RIS channel, theory',...
    'RIS-UE channel, sim.','RIS-UE chhannel, theory');
set(gca,'xminortick','on');
set(gca,'yminortick','on');
grid on;
xlabel('SNR_{pilot}[dB]');
ylabel('NMSE[dB]');