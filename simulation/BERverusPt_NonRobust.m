clc;clear;close all;
cd ..
rng(0);
addpath(genpath(pwd));
NMSE = 0.04; %choose from {0, 0.01, 0.04, 0.1}
Par = load_parameters();
load(['./Data/NMSE=',num2str(NMSE),'ChannelHat.mat'],'Hhat_sd_all','Hhat_sr_all','Hhat_rd_all',...
    'Sigma_sd','Sigma_sr','Sigma_rd',...
    'Psi_sd','Psi_sr','Psi_rd');
load(['./Result/NMSE=',num2str(NMSE),'PerfectCSI_TransceiverResult.mat'],...
    'max_mse_all','B_all','Phi_all','R_all');
Sigma_sd_sqrt = Sigma_sd^(1/2);Psi_sd_sqrt = Psi_sd^(1/2);
Sigma_sr_sqrt = Sigma_sr^(1/2);Psi_sr_sqrt = Psi_sr^(1/2);
Sigma_rd_sqrt = Sigma_rd^(1/2);Psi_rd_sqrt = Psi_rd^(1/2);
%% setting number of antennas
Nt = Par.Nt; 
M = Par.M;
Nr = Par.Nr;
%% number of data streams
D = Par.D;
%% noise power
sigma_q= Par.sigma_q;
%% Pt_dBm_range
Pt_dBm_range = Par.Pt_dBm_range;
%% Monte Carlo
number = Par.DataTransmitionNumber;
par.trials = 1e5; % number of Monte-Carlo trials (transmissions)
%% set up Gray-mapped constellation alphabet
par.mod = 'QPSK'; % modulation type: 'BPSK','QPSK','8PSK','16QAM','64QAM'
par.runId = 0; % simulation ID (used to reproduce results)
% set up Gray-mapped constellation alphabet
switch (par.mod)
    case 'BPSK'
        par.symbols = [ -1 1 ];
    case 'QPSK'
        par.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
    case '8PSK'
        par.symbols = [...
                    exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
                    exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
                    exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
                    exp(1i*2*pi/8*4), exp(1i*2*pi/8*5)];
    case '16QAM'
        par.symbols = [...
                    -3-3i,-3-1i,-3+3i,-3+1i, ...
                    -1-3i,-1-1i,-1+3i,-1+1i, ...
                    +3-3i,+3-1i,+3+3i,+3+1i, ...
                    +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
        par.symbols = [...
                    -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
                    -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
                    -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
                    -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
                    +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
                    +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
                    +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
                    +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    otherwise
        error('constellation not supported!')
end
% normalize symbol energy
par.symbols = par.symbols/sqrt(sum(abs(par.symbols).^2)/length(par.symbols));

% precompute bit labels
par.card = length(par.symbols); % cardinality
par.bps = log2(par.card); % number of bits per symbol
par.bits = de2bi(0:par.card-1,par.bps,'left-msb'); % symbols-to-bits
% initialize result arrays
res.SER = zeros(1,length(Pt_dBm_range));
res.BER  = zeros(1,length(Pt_dBm_range));
res.sumMSE  = zeros(1,length(Pt_dBm_range));
res.per_stream_maxMSE = zeros(1,length(Pt_dBm_range));
% -- start simulation

% trials loop
for num=1:number
    num
    Hhat_sd = Hhat_sd_all(:,:,num);
    Hhat_sr = Hhat_sr_all(:,:,num);
    Hhat_rd = Hhat_rd_all(:,:,num);
    % use runId random seed (enables reproducibility)
    rng(par.runId);
    tmp = zeros(D,length(Pt_dBm_range));
    for tt=1:par.trials

        % generate random bit stream
        b = randi([0 1],D,par.bps);

        % generate transmit symbols
        idx = bi2de(b,'left-msb')+1;
        s = par.symbols(idx).';

        % generate noise vector
        n = sqrt(0.5)*(randn(Nr,1)+1i*randn(Nr,1));
        %  generate total channel
        if NMSE == 0
            Hsd = Hhat_sd;
            Hsr = Hhat_sr;
            Hrd = Hhat_rd;
        else
            Hsd = Hhat_sd+Psi_sd_sqrt*(sqrt(1/2).*(randn(Nr,Nt)+1j*randn(Nr,Nt)))*Sigma_sd_sqrt;
            Hsr = Hhat_sr+Psi_sr_sqrt*(sqrt(1/2).*(randn(M,Nt)+1j*randn(M,Nt)))*Sigma_sr_sqrt;
            Hrd = Hhat_rd+Psi_rd_sqrt*(sqrt(1/2).*(randn(Nr,M)+1j*randn(Nr,M)))*Sigma_rd_sqrt;
        end         
        % SNR loop
        for i = 1:length(Pt_dBm_range)
            B = B_all(:,:,i,num);
            Phi = Phi_all(:,:,i,num);
            R = R_all(:,:,i,num);
            H = Hrd*Phi*Hsr+Hsd;
            % transmit data over noisy channel
            y = H*B*s + sqrt(sigma_q)*n;
            % scale received signal at the UEs (not needed for PSK)
            shat = R'*y; 

            % UE-side nearest-neighbor detection
            [~,idxhat] = min(abs(shat*ones(1,length(par.symbols))-ones(D,1)*par.symbols).^2,[],2); 
            bhat = par.bits(idxhat,:);

            % -- compute error metrics
            err = (idx~=idxhat); % check for symbol errors
            res.SER(i) = res.SER(i) + sum(err)/D; % symbol error rate
            res.BER(i) = res.BER(i) + sum(sum(b~=bhat))/(D*par.bps); % bit error rate
            res.sumMSE(i) = res.sumMSE(i) + sum(abs(shat-s).^2);
            tmp(:,i) = tmp(:,i) + abs(shat-s).^2;
        end
    end
    res.per_stream_maxMSE =res.per_stream_maxMSE+max(tmp/par.trials,[],1);
end
% normalize results
res.SER = res.SER/(par.trials*number);
res.BER = res.BER/(par.trials*number);
res.sumMSE = res.sumMSE/(par.trials*number);
res.per_stream_maxMSE = res.per_stream_maxMSE/number;
save(['./Result/NMSE=',num2str(NMSE),'NonRobustTransceiverSimResult.mat'],'res');
cd simulation