clear
close all
%% Load data
load Tx_sinusoids % = sines(Nf, N_samples), cosines(Nf, N_samples)
load Rx_sines_cable
load Rx_cosines_cable
%% Extract some parameters
N_RxSines = length(Rx_sines); % Number of received sine samples
N_RxCos = length(Rx_cosines); % Number of received cosine samples
Nf = length(sines(:,1)); % Number of sinusoids

%% Try to separate the sinusoids
Rx_sines_mat = zeros(Nf,45000);
skip = 44620;
start = 4.27e4;

Rx_sines_mat(1,1:start) = Rx_sines(1:start);
sample_idx = 1;
for i = 1:1:Nf
   
   
   Rx_sines_mat(i,1:skip) = Rx_sines(idx:idx+skip-1);
    
   sample_idx = sample_idx = sample_idx+1;
end


