% clear
% close all

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
load Tx_sines_5runs_speaker
Tx_sines_5runs_speaker = Tx_sines_5runs_speaker(179010 - dd*Fs: end); % Eyeballing
%181213

load Rx_sines_5runs_speaker; 
Rx_sines_5runs_speaker = Rx_sines_5runs_speaker(136315-dd*Fs:end); % Eyeballing
%136415
load timestamps_5runs



%% Separate the 5 f#@kers
Tx1 = Tx_sines_5runs_speaker(1:Nf*d*Fs);
Rx1 = Rx_sines_5runs_speaker(1:Nf*d*Fs);

Tx2 = Tx_sines_5runs_speaker(Nf*d*Fs +5*Fs - dd*Fs + 1 + 2009        : 2*Nf*d*Fs +5*Fs);
Rx2 = Rx_sines_5runs_speaker(Nf*d*Fs +5*Fs - dd*Fs + 1 + 2009+290    : 2*Nf*d*Fs +5*Fs+290);

Tx3 = Tx_sines_5runs_speaker(2*Nf*d*Fs +10*Fs - dd*Fs + 1 + 1815     : 3*Nf*d*Fs + 10*Fs);
Rx3 = Rx_sines_5runs_speaker(2*Nf*d*Fs +10*Fs - dd*Fs + 1 + 1815+575 : 3*Nf*d*Fs + 10*Fs+575);

Tx4 = Tx_sines_5runs_speaker(3*Nf*d*Fs +15*Fs - dd*Fs + 1 + 1620     : 4*Nf*d*Fs + 15*Fs);
Rx4 = Rx_sines_5runs_speaker(3*Nf*d*Fs +15*Fs - dd*Fs + 1 + 1620+865 : 4*Nf*d*Fs + 15*Fs+865);

Tx5 = Tx_sines_5runs_speaker(4*Nf*d*Fs +20*Fs - dd*Fs + 1 + 1426      : end);
Rx5 = Rx_sines_5runs_speaker(4*Nf*d*Fs +20*Fs - dd*Fs + 1 + 1426+1155 : end);

% % figure
% % plot(Tx5); hold on;
% % plot(Rx5.*50);
% % return
%% Separate sines from all of them
Nfft = 2^15;


Tx_mat1 = zeros(Nf, Fs*d); Rx_mat1 = zeros(Nf, Fs*d);
Tx_mat2 = zeros(Nf, Fs*d); Rx_mat2 = zeros(Nf, Fs*d);
Tx_mat3 = zeros(Nf, Fs*d); Rx_mat3 = zeros(Nf, Fs*d);
Tx_mat4 = zeros(Nf, Fs*d); Rx_mat4 = zeros(Nf, Fs*d);
Tx_mat5 = zeros(Nf, Fs*d); Rx_mat5 = zeros(Nf, Fs*d);


Tx1_F = zeros(Nf, Nfft); Rx1_F = zeros(Nf, Nfft);
Tx2_F = zeros(Nf, Nfft); Rx2_F = zeros(Nf, Nfft);
Tx3_F = zeros(Nf, Nfft); Rx3_F = zeros(Nf, Nfft);
Tx4_F = zeros(Nf, Nfft); Rx4_F = zeros(Nf, Nfft);
Tx5_F = zeros(Nf, Nfft); Rx5_F = zeros(Nf, Nfft);

mag_response1 = zeros(1,Nf); phase_response1 = zeros(1,Nf);
mag_response2 = zeros(1,Nf); phase_response2 = zeros(1,Nf);
mag_response3 = zeros(1,Nf); phase_response3 = zeros(1,Nf);
mag_response4 = zeros(1,Nf); phase_response4 = zeros(1,Nf);
mag_response5 = zeros(1,Nf); phase_response5 = zeros(1,Nf);


win_idx = 4000:12000; % Time domain window for DFT (Eyeballing)

for i = 1:1:Nf
    % Separate them
    idx = (i-1)*Fs*d + 1;
    Tx_mat1(i,:) = Tx1(idx:idx + Fs*d -1);  Rx_mat1(i,:) = Rx1(idx:idx + Fs*d -1);
    Tx_mat2(i,:) = Tx2(idx:idx + Fs*d -1);  Rx_mat2(i,:) = Rx2(idx:idx + Fs*d -1);
    Tx_mat3(i,:) = Tx3(idx:idx + Fs*d -1);  Rx_mat3(i,:) = Rx3(idx:idx + Fs*d -1);
    Tx_mat4(i,:) = Tx4(idx:idx + Fs*d -1);  Rx_mat4(i,:) = Rx4(idx:idx + Fs*d -1);
    Tx_mat5(i,:) = Tx5(idx:idx + Fs*d -1);  Rx_mat5(i,:) = Rx5(idx:idx + Fs*d -1);
    
    % Get them in freq
    [Tx1_F(i,:) , w] = freqz(Tx_mat1(i,win_idx),1,2^15,Fs); [Rx1_F(i,:) , w] = freqz(Rx_mat1(i,win_idx),1,2^15,Fs);
    [Tx2_F(i,:) , w] = freqz(Tx_mat2(i,win_idx),1,2^15,Fs); [Rx2_F(i,:) , w] = freqz(Rx_mat2(i,win_idx),1,2^15,Fs);
    [Tx3_F(i,:) , w] = freqz(Tx_mat3(i,win_idx),1,2^15,Fs); [Rx3_F(i,:) , w] = freqz(Rx_mat3(i,win_idx),1,2^15,Fs);
    [Tx4_F(i,:) , w] = freqz(Tx_mat4(i,win_idx),1,2^15,Fs); [Rx4_F(i,:) , w] = freqz(Rx_mat4(i,win_idx),1,2^15,Fs);
    [Tx5_F(i,:) , w] = freqz(Tx_mat5(i,win_idx),1,2^15,Fs); [Rx5_F(i,:) , w] = freqz(Rx_mat5(i,win_idx),1,2^15,Fs);
    
    % Find the closest value of freq_arr(i) inside w
    [w_idx, w_idx] = min(abs(w - freq_arr(i))); 
 
    %% Find the max value 'in the vicinity' (+-delta) of the corresponding frequency
    % Make an index array of length +- delta
    delta = 50;
    f_idx = w_idx-delta:w_idx+delta;
    if w_idx+delta >= Nfft
        f_idx = w_idx-delta:1:Nfft;
    end
    
    % Find index where max happens for the 5 experiments
    index_maxTx1 = find(abs(Tx1_F(i,f_idx)) == max(abs(Tx1_F(i,f_idx)))) + w_idx - delta-1;
    index_maxRx1 = find(abs(Rx1_F(i,f_idx)) == max(abs(Rx1_F(i,f_idx)))) + w_idx - delta-1;
    
    index_maxTx2 = find(abs(Tx2_F(i,f_idx)) == max(abs(Tx2_F(i,f_idx)))) + w_idx - delta-1;
    index_maxRx2 = find(abs(Rx2_F(i,f_idx)) == max(abs(Rx2_F(i,f_idx)))) + w_idx - delta-1;
    
    index_maxTx3 = find(abs(Tx3_F(i,f_idx)) == max(abs(Tx3_F(i,f_idx)))) + w_idx - delta-1;
    index_maxRx3 = find(abs(Rx3_F(i,f_idx)) == max(abs(Rx3_F(i,f_idx)))) + w_idx - delta-1;
    
    index_maxTx4 = find(abs(Tx4_F(i,f_idx)) == max(abs(Tx4_F(i,f_idx)))) + w_idx - delta-1;
    index_maxRx4 = find(abs(Rx4_F(i,f_idx)) == max(abs(Rx4_F(i,f_idx)))) + w_idx - delta-1;
    
    index_maxTx5 = find(abs(Tx5_F(i,f_idx)) == max(abs(Tx5_F(i,f_idx)))) + w_idx - delta-1;
    index_maxRx5 = find(abs(Rx5_F(i,f_idx)) == max(abs(Rx5_F(i,f_idx)))) + w_idx - delta-1;
    
    %% Get mag and phase response for the 5 experiments
    mag_response1(i) = db(Rx1_F(i,index_maxRx1)) - db(Tx1_F(i,index_maxTx1)); % dB
    phase_response1(i) = angle(Rx1_F(i,index_maxRx1)/Tx1_F(i,index_maxTx1)).*(180/(pi)); % degrees
    
    mag_response2(i) = db(Rx2_F(i,index_maxRx2)) - db(Tx2_F(i,index_maxTx2));
    phase_response2(i) = angle(Rx2_F(i,index_maxRx2)/Tx2_F(i,index_maxTx2)).*(180/(pi));
    
    mag_response3(i) = db(Rx3_F(i,index_maxRx3)) - db(Tx3_F(i,index_maxTx3));
    phase_response3(i) = angle(Rx3_F(i,index_maxRx3)/Tx3_F(i,index_maxTx3)).*(180/(pi));
    
    mag_response4(i) = db(Rx4_F(i,index_maxRx4)) - db(Tx4_F(i,index_maxTx4));
    phase_response4(i) = angle(Rx4_F(i,index_maxRx4)/Tx4_F(i,index_maxTx4)).*(180/(pi));
    
    mag_response5(i) = db(Rx5_F(i,index_maxRx5)) - db(Tx5_F(i,index_maxTx5));
    phase_response5(i) = angle(Rx5_F(i,index_maxRx5)/Tx5_F(i,index_maxTx5)).*(180/(pi));
    
end
%% Look at the separated Tx and Rx sines in frequency

figure;
for i = 1:1:Nf
   plot(w,db(Tx1_F(i,:))); hold on; 
end
title('Send (Tx) sinusoids');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold off;

% figure;
% for i = 1:1:Nf
%    plot(w,db(Rx1_F(i,:))); hold on; 
% end
% title('Received (Rx) sinusoids');
% xlabel('Frequency [Hz]');
% ylabel('Magnitude [dB]');
% hold off;

%% Mag and phase of the 5 experiments
figure;
plot(freq_arr/1000,mag_response1); hold on;
plot(freq_arr/1000,mag_response2);hold on;
plot(freq_arr/1000,mag_response3); hold on;
plot(freq_arr/1000,mag_response4); hold on;
plot(freq_arr/1000,mag_response5); 
title('Magnitude response of mic/speaker channel (measured 5 times)');
xlabel('Frequency [kHz]');
ylabel('Magnitude [dB]');
legend('1','2','3','4','5');

%Phase
figure;
plot(freq_arr/1000,unwrap(phase_response1.*pi/180)); hold on;
plot(freq_arr/1000,unwrap(phase_response2.*pi/180)); 
plot(freq_arr/1000,unwrap(phase_response3.*pi/180)); 
plot(freq_arr/1000,unwrap(phase_response4.*pi/180)); 
plot(freq_arr/1000,unwrap(phase_response5.*pi/180)); 
legend('1','2','3','4','5');
title('Phase response');
xlabel('Frequency [kHz]');
ylabel('Phase [Degrees]');