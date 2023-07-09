clc;clear;close all;
rng(0);
warning('off');
addpath(genpath(pwd));
NMSE = 0.04;
Par = load_parameters();
load(['./Data/NMSE=',num2str(NMSE),'ChannelHat_Mmax64.mat'],'Hhat_sd_all','Hhat_sr_all','Hhat_rd_all',...
    'Sigma_sd','Sigma_sr','Sigma_rd',...
    'Psi_sd','Psi_sr','Psi_rd');
Mrange = [4,8,16,32,64];
%% setting number of antennas
Nt = Par.Nt; 
Mmax = Par.Mmax;
Nr = Par.Nr;
%% number of data streams
D = Par.D;
%% parameters of algorithm
iter_max = Par.iter_max;
epsilon = Par.epsilon;
%% noise power
sigma_q= Par.sigma_q;
%% Pt_dBm
Pt_dBm = Par.Pt_dBm; 
Pt = 10^(Pt_dBm/10-3);
%% Monte Carlo
number = Par.DataTransmitionNumber;
%% Initialize
phi_all = exp(1i*2*pi*rand(Mmax,length(Mrange),number));
%% store result
max_mse_all = zeros(length(Mrange),number);
B_all = zeros(Nt,D,length(Mrange),number);
Phi_all = zeros(Mmax,Mmax,length(Mrange),number);
R_all = zeros(Nr,D,length(Mrange),number);
%% channel realizations
for num = 1:number
    Hhat_sd = Hhat_sd_all(:,:,num);
    Hhat_sr = Hhat_sr_all(:,:,num);
    Hhat_rd = Hhat_rd_all(:,:,num);
    %% number of RIS elements
    for i = 1:length(Mrange)
        M = Mrange(i);
        Hhat_sr0 = Hhat_sr(1:M,:);
        Hhat_rd0 = Hhat_rd(:,1:M);
        phi = phi_all(1:M,i,num);
        Sigma_rd0 = Sigma_rd(1:M,1:M);
        Psi_sr0 = Psi_sr(1:M,1:M);
        [max_mse,B,Phi,R] = ManoptSolver(Hhat_sd,Hhat_sr0,Hhat_rd0,...
            Sigma_sd,Sigma_sr,Sigma_rd0,Psi_sd,Psi_sr0,Psi_rd,...
            Pt,sigma_q,Nt,M,Nr,D,phi,iter_max,epsilon);
        max_mse_all(i,num) = max_mse;
        B_all(:,:,i,num) = B;
        Phi_all(1:M,1:M,i,num) = Phi;
        R_all(:,:,i,num) = R;
    end
end
save(['./Result/NMSE=',num2str(NMSE),'MSEverusM_TransceiverResult.mat'],...
    'max_mse_all','B_all','Phi_all','R_all');