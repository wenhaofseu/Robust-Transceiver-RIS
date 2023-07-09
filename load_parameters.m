function Par = load_parameters()
    %% setting number of antennas
    Par.Nt = 16; %BS
    Par.M = 32; %RIS(default)
    Par.Nr = 4; %UE
    Par.Mmax = 128; %RIS(maximum)
    %%  setting correlation coefficient, i.e., \rho
    Par.rho_BS = 0.6; 
    Par.rho_RIS_Tx = 0.8; Par.rho_RIS_Rx = 0.8;
    Par.rho_UE = 0.4;
    %% setting pass loss exponent and reference path loss
    Par.alpha_sd=4.0; Par.PL0_sd = 35; %BS--UE link
    Par.alpha_sr=2.0; Par.PL0_sr = 30; %BS--RIS link
    Par.alpha_rd=2.0; Par.PL0_rd = 30; %RIS--UE link
    %% setting Rician factor
    Par.kappa_sd = 0; %BS--UE link
    Par.kappa_sr = 5; %BS--RIS link
    Par.kappa_rd = 5; %RIS--UE link
    %% setting location 
    Par.BS = [0,0];
    Par.RIS = [50,5];
    Par.UE = [65,0];
    %% Monte Carlo
    Par.ChannelEstimationNumber = 1e4; %Monte Carlo for Channel Estimation
    Par.DataTransmitionNumber = 2e3; %Monte Carlo for Data Transmition
    %% noise power
    Par.sigma_q = 1e-11; %-80dBm
    %% SNR_range for channel estimation
    Par.SNR_range = 0:5:30;
    %% Pt_dBm_range
    Par.Pt_dBm_range = [15:4:23,25:2:35];
    Par.Pt_dBm = 25;
    %% number of data streams
    Par.D = 3;
    %% parameters of algorithm
    Par.iter_max = 1000;
    Par.epsilon = 1e-3;
    %% parameter of GaussianRandomization
    Par.GaussianRandomizationNumber = 1000;
end
