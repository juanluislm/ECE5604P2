clear
close all
%% Parameters
% Sampling frequency
Fs = 44.1e3; 
% Duration of each sinusoid
dd = 0.05;
d = 0.25 + dd; 
df = 0.005;
% Array of frequencies
freq_arr = (df:df:0.5).*Fs; 
Nf = length(freq_arr) % Length of frequency array
% Time domain vector
nT = (0:d*Fs)/Fs;
N_samples = length(nT); % Length of time domain vector


%% Create Session
s = daq.createSession('digilent')
%% Add IO channels
% The Digilent Analog Discovery board analog channels are labeled starting
% from '1':
% Analog Inputs:    '1+'/'1-' or '2+/2-'
% Analog Outputs:   'W1' / 'W2'.
ch = s.addAnalogInputChannel('AD1',1, 'Voltage');
%s.addAnalogOutputChannel('AD1',1, 'Voltage');
%s.Channels

s.Rate = 44.1e3;  % set rate of object to desired rate
s.Channels.Range = [-5 5];
s.DurationInSeconds = Nf*d + 5;

%% 


%% Create the sampled sines and cosines
sines = zeros(Nf, N_samples);
cosines = zeros(Nf, N_samples);
for i = 1:1:Nf
    sines(i,:) = sin(2*pi*freq_arr(i).*nT);
    cosines(i,:) = cos(2*pi*freq_arr(i).*nT);
end
[n,m] = size(sines); 
sines(:,m-(dd*Fs)+1:m) = zeros(Nf,dd*Fs); cosines(:,m-(dd*Fs)+1:m) = zeros(Nf,dd*Fs);
sines_col = sines'; cosines_col = cosines';


sines = reshape(sines_col,1,numel(sines_col));
cosines = reshape(cosines_col,1,numel(cosines_col));

sines_Z = [zeros(1,5*Fs) , sines];
cosines_Z = [zeros(1,5*Fs) , cosines];

%% Send
sound(sines_Z,Fs);
[DigilentData, timestamps, triggerTime] = startForeground(s);
%pause(d*2)
Tx_sines_after_soundCard_headSet = DigilentData;
Tx_sines_before_soundCard_headSet = sines_Z;

save('Tx_sines_after_soundCard_headSet.mat','Tx_sines_after_soundCard_headSet');
save('Tx_sines_before_soundCard_headSet.mat','Tx_sines_before_soundCard_headSet');

% Tx_cosines_after_soundCard_headSet = DigilentData;
% Tx_cosines_before_soundCard_headSet = cosines_Z;
% 
% save('Tx_cosines_after_soundCard_headSet.mat','Tx_cosines_after_soundCard_headSet');
% save('Tx_cosines_before_soundCard_headSet.mat','Tx_cosines_before_soundCard_headSet');

% 
% save('timestamps.mat','timestamps')
% clear sound
