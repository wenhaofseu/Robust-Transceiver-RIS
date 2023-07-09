clc;clear;close all;
rng(0);
warning('off');
addpath(genpath(pwd));
NMSE = 0.04; %choose from {0, 0.01, 0.04, 0.1}
Par = load_parameters();
load(['./Data/NMSE=',num2str(NMSE),'ChannelHat.mat'],'Hhat_sd_all','Hhat_sr_all','Hhat_rd_all',...
    'Sigma_sd','Sigma_sr','Sigma_rd',...
    'Psi_sd','Psi_sr','Psi_rd');
%% setting number of antennas
Nt = Par.Nt; 
M = Par.M;
Nr = Par.Nr;
%% number of data streams
D = Par.D;
%% parameters of algorithm
iter_max = Par.iter_max;
epsilon = Par.epsilon;
%% noise power
sigma_q= Par.sigma_q;
%% Pt_dBm_range
Pt_dBm_range = Par.Pt_dBm_range;
%% Monte Carlo
number = Par.DataTransmitionNumber;
phi_all = exp(1i*2*pi*rand(M,1,length(Pt_dBm_range),number));
%% store result
max_mse_all = zeros(length(Pt_dBm_range),number);
B_all = zeros(Nt,D,length(Pt_dBm_range),number);
Phi_all = zeros(M,M,length(Pt_dBm_range),number);
R_all = zeros(Nr,D,length(Pt_dBm_range),number);
for num = 1:number
    num
    Hhat_sd = Hhat_sd_all(:,:,num);
    Hhat_sr = Hhat_sr_all(:,:,num);
    Hhat_rd = Hhat_rd_all(:,:,num);
    for i = 1:length(Pt_dBm_range)
        i
        Pt_dBm = Pt_dBm_range(i);
        Pt = 10^(Pt_dBm/10-3);
        phi = phi_all(:,:,i,num);
        [max_mse,B,Phi,R] = ManoptSolver(Hhat_sd,Hhat_sr,Hhat_rd,...
            Sigma_sd,Sigma_sr,Sigma_rd,Psi_sd,Psi_sr,Psi_rd,...
            Pt,sigma_q,Nt,M,Nr,D,phi,iter_max,epsilon);
        max_mse_all(i,num) = max_mse;
        B_all(:,:,i,num) = B;
        Phi_all(:,:,i,num) = Phi;
        R_all(:,:,i,num) = R;
    end
end
save(['./Result/NMSE=',num2str(NMSE),'TransceiverResult.mat'],...
    'max_mse_all','B_all','Phi_all','R_all');