clear
close all

%% Parameters
% Sampling frequency
Fs = 44.1e3; 
% Duration of each sinusoid
dd = 0.05;
d = 0.25 + dd; 
% Array of frequencies
df = 0.005;
freq_arr = (df:df:0.5).*Fs; 
Nf = length(freq_arr); % Length of frequency array
% Time domain vector
nT = (0:d*Fs)/Fs;
N_samples = length(nT); % Length of time domain vector
%% Load data
load Tx_sines_after_soundCard 
Tx_sines_after_soundCard = Tx_sines_after_soundCard(181850 - dd*Fs: end); % Eyeballing
load Rx_sines_tapped_measurement; 
Rx_sines_tapped_measurement = Rx_sines_tapped_measurement(121465-dd*Fs:end); % Eyeballing
%121380

load Tx_sines_after_soundCard_headSet 
load Rx_sines_headset %Rx_sines_headset = Rx_sines_headset(140700:end); % Eyeballing

Tx_sines_after_soundCard_headSet = Tx_sines_after_soundCard_headSet(181850 - dd*Fs: end); % Eyeballing
Rx_sines_headset = Rx_sines_headset(127508-dd*Fs:end);%116750
%for idx=1 -> 127508
%for idx=2 -> 127563
%for idx=3 -> 127578
%for idx=2 -> 127563

%% Extract some parameters

Rx_sines_mat_speaker = zeros(Nf, Fs*d);
Tx_sines_mat_speaker = zeros(Nf, Fs*d);

Nfft = 2^15;
Tx_F_speaker = zeros(Nf, Nfft);
Rx_F_speaker = zeros(Nf, Nfft);
mag_response_speaker = zeros(1,Nf);
phase_response_speaker = zeros(1,Nf);
win_idx = 4000:12000; % Eyeballing

Rx_sines_mat_headset = zeros(Nf, Fs*d);
Tx_sines_mat_headset = zeros(Nf, Fs*d);

Nfft = 2^15;
Tx_F_headset = zeros(Nf, Nfft);
Rx_F_headset = zeros(Nf, Nfft);
mag_response_headset = zeros(1,Nf);
phase_response_headset = zeros(1,Nf);

for i = 1:1:Nf
    % Separate them
%     idx = (i-1)*Fs*d + 1;
%     Rx_sines_mat_speaker(i,:) = Rx_sines_tapped_measurement(idx:idx + Fs*d -1);
%     Tx_sines_mat_speaker(i,:) = Tx_sines_after_soundCard(idx:idx + Fs*d -1);
%     
%     % Get them in freq
%     [Tx_F_speaker(i,:) , w] = freqz(Tx_sines_mat_speaker(i,win_idx),1,2^15,Fs);
%     [Rx_F_speaker(i,:) , w] = freqz(Rx_sines_mat_speaker(i,win_idx),1,2^15,Fs);
%     
%     % Find the closest value of freq_arr(i) inside w
%     [w_idx, w_idx] = min(abs(w - freq_arr(i))); 
%  
%     % Find the max value 'in the vicinity' (+-delta) of the corresponding frequency
%     delta = 50;
%     f_idx = w_idx-delta:w_idx+delta;
%     if w_idx+delta >= Nfft
%         f_idx = w_idx-delta:1:Nfft;
%     end
%     
%     index_maxTx = find(abs(Tx_F_speaker(i,f_idx)) == max(abs(Tx_F_speaker(i,f_idx)))) + w_idx - delta-1;
%     index_maxRx = find(abs(Rx_F_speaker(i,f_idx)) == max(abs(Rx_F_speaker(i,f_idx)))) + w_idx - delta-1;
%     
%     % Get mag and phase response
%     mag_response_speaker(i) = db(Rx_F_speaker(i,index_maxRx)) - db(Tx_F_speaker(i,index_maxTx));
%     phase_response_speaker(i) = angle(Rx_F_speaker(i,index_maxRx)/Tx_F_speaker(i,index_maxTx)).*(180/(pi));
    
    %% Headset
    idx = (i-1)*Fs*d + 1;
    Rx_sines_mat_headset(i,:) = Rx_sines_headset(idx:idx + Fs*d -1);
    Tx_sines_mat_headset(i,:) = Tx_sines_after_soundCard_headSet(idx:idx + Fs*d -1);
    
    % Get them in freq
    [Tx_F_headset(i,:) , w] = freqz(Tx_sines_mat_headset(i,win_idx),1,2^15,Fs);
    [Rx_F_headset(i,:) , w] = freqz(Rx_sines_mat_headset(i,win_idx),1,2^15,Fs);
    
    % Find the closest value of freq_arr(i) inside w
    [w_idx, w_idx] = min(abs(w - freq_arr(i))); 
 
    % Find the max value 'in the vicinity' (+-delta) of the corresponding frequency
    delta = 100;
    f_idx = w_idx-delta:w_idx+delta;
    if w_idx+delta >= Nfft
        f_idx = w_idx-delta:1:Nfft;
    end
    
    index_maxTx = find(abs(Tx_F_headset(i,f_idx)) == max(abs(Tx_F_headset(i,f_idx)))) + w_idx - delta-1;
    index_maxRx = find(abs(Rx_F_headset(i,f_idx)) == max(abs(Rx_F_headset(i,f_idx)))) + w_idx - delta-1;
    
    % Get mag and phase response
    mag_response_headset(i) = db(Rx_F_headset(i,index_maxRx)) - db(Tx_F_headset(i,index_maxTx));
    phase_response_headset(i) = angle(Rx_F_headset(i,index_maxRx)/Tx_F_headset(i,index_maxTx)).*(180/(pi));
    
end

plot_idx = 3;
scale = max(abs(Tx_sines_mat_headset(plot_idx,:)))/max(abs(Rx_sines_mat_headset(plot_idx,:)));
% Time
figure; 
plot(Tx_sines_mat_headset(plot_idx,:)); hold on; plot(Rx_sines_mat_headset(plot_idx,:).*scale,'r');

return

% plot_idx = 70;
% scale = max(abs(Tx_sines_mat_speaker(plot_idx,:)))/max(abs(Rx_sines_mat_speaker(plot_idx,:)));
% % % Time
% % figure; 
% % plot(Tx_sines_mat_speaker(plot_idx,:)); hold on; plot(Rx_sines_mat_speaker(plot_idx,:).*scale);
% % 
% % % Freq
% % figure; 
% % plot(w,db(Tx_F_speaker(plot_idx,:))); hold on; plot(w, db(Rx_F_speaker(plot_idx,:).*scale));
% 
% 
% figure;
% for i = 1:1:Nf
%    plot(w,db(Tx_F_speaker(i,:))); hold on; 
% end
% title('Send (Tx) sinusoids');
% xlabel('Frequency [Hz]');
% ylabel('Magnitude [dB]');
% hold off;
% 
% figure;
% for i = 1:1:Nf
%    plot(w,db(Rx_F_speaker(i,:))); hold on; 
% end
% title('Received (Rx) sinusoids');
% xlabel('Frequency [Hz]');
% ylabel('Magnitude [dB]');
% hold off;
% 
% 
% figure;
% plot(freq_arr,mag_response_speaker);
% title('Magnitude response');
% xlabel('Frequency [Hz]');
% ylabel('Magnitude [dB]');
% 
% figure;
% plot(freq_arr,phase_response_speaker);
% title('Phase response');
% xlabel('Frequency [Hz]');
% ylabel('Phase [Degrees]');