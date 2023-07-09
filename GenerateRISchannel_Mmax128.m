clc;clear;close all;
%% setting seed for reproduction
rng(0);
%% the function 'load_parameters' needs to be called
Par = load_parameters();
%% setting number of antennas
Nt = Par.Nt; 
M = Par.Mmax;
Nr = Par.Nr; 
%%  setting correlation coefficient
rho_BS = Par.rho_BS; 
rho_RIS_Tx = Par.rho_RIS_Tx; rho_RIS_Rx = Par.rho_RIS_Rx;
rho_UE = Par.rho_UE;
%% setting path loss exponent and reference path loss
alpha_sd = Par.alpha_sd; PL0_sd = Par.PL0_sd;
alpha_sr = Par.alpha_sr; PL0_sr = Par.PL0_sr;
alpha_rd = Par.alpha_rd; PL0_rd = Par.PL0_rd;
%% setting Rician factor
kappa_sd = Par.kappa_sd;
kappa_sr = Par.kappa_sr;
kappa_rd = Par.kappa_rd;
%% setting location 
BS = Par.BS;
RIS = Par.RIS;
UE = Par.UE;
%% Monte Carlo
number = Par.ChannelEstimationNumber;
%%  Generate covariance matrix according to correlation coefficient
R_Tx_sd = covariance_matrix(rho_BS,Nt); R_Rx_sd = covariance_matrix(rho_UE,Nr);
R_Tx_sr = covariance_matrix(rho_BS,Nt); R_Rx_sr = covariance_matrix(rho_RIS_Rx,M);
R_Tx_rd = covariance_matrix(rho_RIS_Tx,M); R_Rx_rd = covariance_matrix(rho_UE,Nr);
%% Precompute the square root of covariance matrix
R_Tx_sd_sqrt = R_Tx_sd^(1/2); R_Rx_sd_sqrt = R_Rx_sd^(1/2);
R_Tx_sr_sqrt = R_Tx_sr^(1/2); R_Rx_sr_sqrt = R_Rx_sr^(1/2);
R_Tx_rd_sqrt = R_Tx_rd^(1/2); R_Rx_rd_sqrt = R_Rx_rd^(1/2);
%% Calculate distance [meter]
d_sd = norm(BS-UE);
d_sr = norm(BS-RIS);
d_rd = norm(RIS-UE);
%% path_loss 
path_loss = @(d,alpha,PL0) -PL0-10*alpha*log10(d); %path loss formula
PL_sd = path_loss(d_sd,alpha_sd,PL0_sd);
PL_sr = path_loss(d_sr,alpha_sr,PL0_sr);
PL_rd = path_loss(d_rd,alpha_rd,PL0_rd);
pl_sd = 10^(PL_sd/10);
pl_sr = 10^(PL_sr/10);
pl_rd = 10^(PL_rd/10);
%% Generate AOD and AOA
AOD_sd = pi*rand-pi/2; AOA_sd = pi*rand-pi/2;
AOD_sr = pi*rand-pi/2; AOA_sr = pi*rand-pi/2;
AOD_rd = pi*rand-pi/2; AOA_rd = pi*rand-pi/2;
%% calculate LOS component
Hbar_sd = ULA_SteeringVector(AOA_sd,Nr)*ULA_SteeringVector(AOD_sd,Nt)';
Hbar_sr = ULA_SteeringVector(AOA_sr,M)*ULA_SteeringVector(AOD_sr,Nt)';
Hbar_rd = ULA_SteeringVector(AOA_rd,Nr)*ULA_SteeringVector(AOD_rd,M)';
%% calculate radio of LOS and NLOS
kappa2_sd = 1/(1+kappa_sd); kappa1_sd = 1-kappa2_sd;
kappa1_sd = sqrt(kappa1_sd); kappa2_sd = sqrt(kappa2_sd);
kappa2_sr = 1/(1+kappa_sr); kappa1_sr = 1-kappa2_sr;
kappa1_sr = sqrt(kappa1_sr); kappa2_sr = sqrt(kappa2_sr);
kappa2_rd = 1/(1+kappa_rd); kappa1_rd = 1-kappa2_rd;
kappa1_rd = sqrt(kappa1_rd); kappa2_rd = sqrt(kappa2_rd);
%% store matrix
Hsd_all = zeros(Nr,Nt,number);
Hsr_all = zeros(M,Nt,number);
Hrd_all = zeros(Nr,M,number);
%% setting seed for reproduction
rng(0);
%% channel realization
for num = 1:number
	%% generate NLOS
    Htilde_sd = sqrt(1/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
    Htilde_sr = sqrt(1/2)*(randn(M,Nt)+1i*randn(M,Nt));
    Htilde_rd = sqrt(1/2)*(randn(Nr,M)+1i*randn(Nr,M));
	%% combine LOS and NLOS 
    Hsd = sqrt(pl_sd)*(kappa1_sd*Hbar_sd+kappa2_sd*R_Rx_sd_sqrt*Htilde_sd*R_Tx_sd_sqrt);
    Hsr = sqrt(pl_sr)*(kappa1_sr*Hbar_sr+kappa2_sr*R_Rx_sr_sqrt*Htilde_sr*R_Tx_sr_sqrt);
    Hrd = sqrt(pl_rd)*(kappa1_rd*Hbar_rd+kappa2_rd*R_Rx_rd_sqrt*Htilde_rd*R_Tx_rd_sqrt);
    Hsd_all(:,:,num) = Hsd;
    Hsr_all(:,:,num) = Hsr;
    Hrd_all(:,:,num) = Hrd;
end
save('RISchannel_Mmax128.mat','Hsd_all','Hsr_all','Hrd_all',...
    'Hbar_sd','Hbar_sr','Hbar_rd',...
    'R_Tx_sd','R_Rx_sd','R_Tx_sr','R_Rx_sr','R_Tx_rd','R_Rx_rd',...
    'pl_sd','pl_sr','pl_rd');

function R = covariance_matrix(rho,N)
    R = zeros(N);
    index = 1:N;
    for i = index
        for j = index
            R(i,j) = rho^abs(i-j);
        end
    end
end

function h = ULA_SteeringVector(phi,N)
    h = exp(1i*pi*sin(phi)*(0:N-1)'); %element space equals lambda/2
end