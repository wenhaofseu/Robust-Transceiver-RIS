clc;clear;close all;
NMSE = 0.04; %choose from {0, 0.01, 0.04, 0.1}
Par = load_parameters();
load('RISchannel.mat','Hsd_all','Hsr_all','Hrd_all',...
    'Hbar_sd','Hbar_sr','Hbar_rd',...
    'R_Tx_sd','R_Rx_sd','R_Tx_sr','R_Rx_sr','R_Tx_rd','R_Rx_rd',...
    'pl_sd','pl_sr','pl_rd');
%% setting number of antennas
Nt = Par.Nt; 
M = Par.M;
Nr = Par.Nr;
%% Precompute Inv
Inv_R_Rx_sd = inv(R_Rx_sd);
Inv_R_Rx_sr = inv(R_Rx_sr);
Inv_R_Rx_rd = inv(R_Rx_rd);
%% Precompute square root
R_Tx_sd_sqrt = R_Tx_sd^(1/2); R_Rx_sd_sqrt = R_Rx_sd^(1/2);
R_Tx_sr_sqrt = R_Tx_sr^(1/2); R_Rx_sr_sqrt = R_Rx_sr^(1/2);
R_Tx_rd_sqrt = R_Tx_rd^(1/2); R_Rx_rd_sqrt = R_Rx_rd^(1/2);
%% Get \Sigma
Sigma_sd = R_Tx_sd;
Sigma_sr = R_Tx_sr;
Sigma_rd = R_Tx_rd;
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
%% Monte Carlo
number = Par.DataTransmitionNumber;
%% Calculate the pilot signal-to-noise ratio based on NMSE.
NMSE_sd = NMSE;
func = @(sigma_0_q_sd) sigma_0_q_sd*trace(inv(eye(Nr)+sigma_0_q_sd.*Inv_R_Rx_sd))/(1+kappa_sd)/Nr-NMSE_sd;
[sigma_0_q_sd,~] = fsolve(func,0);
NMSE_sr = NMSE;
func = @(sigma_0_q_sr) sigma_0_q_sr*trace(inv(eye(M)+sigma_0_q_sr.*Inv_R_Rx_sr))/(1+kappa_sr)/Nr-NMSE_sr;
[sigma_0_q_sr,~] = fsolve(func,0);
NMSE_rd = NMSE;
func = @(sigma_0_q_rd) sigma_0_q_rd*trace(inv(eye(Nr)+sigma_0_q_rd.*Inv_R_Rx_rd))/(1+kappa_rd)/Nr-NMSE_rd;
[sigma_0_q_rd,~] = fsolve(func,0);
sigma_ce_q_sd = sigma_0_q_sd*pl_sd/(1+kappa_sd);
sigma_ce_q_sr = sigma_0_q_sr*pl_sr/(1+kappa_sr);
sigma_ce_q_rd = sigma_0_q_rd*pl_rd/(1+kappa_rd);
Psi_sd = sigma_ce_q_sd.*inv(eye(Nr)+sigma_0_q_sd.*Inv_R_Rx_sd);
Psi_sr = sigma_ce_q_sr.*inv(eye(M)+sigma_0_q_sr.*Inv_R_Rx_sr);
Psi_rd = sigma_ce_q_rd.*inv(eye(Nr)+sigma_0_q_rd.*Inv_R_Rx_rd);
%% calculate Temp
Temp_sd = (eye(Nr)+sigma_0_q_sd.*Inv_R_Rx_sd)^(-1/2);
Temp_sr = (eye(M)+sigma_0_q_sr.*Inv_R_Rx_sr)^(-1/2);
Temp_rd = (eye(Nr)+sigma_0_q_rd.*Inv_R_Rx_rd)^(-1/2);
%% store matrix
Hhat_sd_all = zeros(Nr,Nt,number);
Hhat_sr_all = zeros(M,Nt,number);
Hhat_rd_all = zeros(Nr,M,number);
%% setting seed for reproduction
rng(0);
for num = 1:number
    Htilde_sd = sqrt(1/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
    Htilde_sr = sqrt(1/2)*(randn(M,Nt)+1i*randn(M,Nt));
    Htilde_rd = sqrt(1/2)*(randn(Nr,M)+1i*randn(Nr,M));
    Hhat_sd = sqrt(pl_sd)*(kappa1_sd*Hbar_sd+kappa2_sd*R_Rx_sd_sqrt*Temp_sd*Htilde_sd*R_Tx_sd_sqrt);
    Hhat_sr = sqrt(pl_sr)*(kappa1_sr*Hbar_sr+kappa2_sr*R_Rx_sr_sqrt*Temp_sr*Htilde_sr*R_Tx_sr_sqrt);
    Hhat_rd = sqrt(pl_rd)*(kappa1_rd*Hbar_rd+kappa2_rd*R_Rx_rd_sqrt*Temp_rd*Htilde_rd*R_Tx_rd_sqrt);
    Hhat_sd_all(:,:,num) = Hhat_sd;
    Hhat_sr_all(:,:,num) = Hhat_sr;
    Hhat_rd_all(:,:,num) = Hhat_rd;
end
save(['./Data/NMSE=',num2str(NMSE),'ChannelHat.mat'],...
	'Hhat_sd_all','Hhat_sr_all','Hhat_rd_all',...
    'Sigma_sd','Sigma_sr','Sigma_rd',...
    'Psi_sd','Psi_sr','Psi_rd');