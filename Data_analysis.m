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
thresh = 0.5e-3;
buff = ones(1,500);
Rx_sines_mat = zeros(Nf,round(N_RxSines/Nf)+100e3);
% For the sines:
freq_idx_counter = 1;
sample_idx_counter = 1;
for i = 1:1:N_RxSines
    Rx_sines_mat(freq_idx_counter,sample_idx_counter) = Rx_sines(i);
    sample_idx_counter = sample_idx_counter+1;
    
    
    buff = [buff(2:end),Rx_sines(i)];
    if mean(abs(buff)) < 5e-4
        freq_idx_counter = freq_idx_counter+1;
        buff = ones(1,1500);
    end
    
end
 
