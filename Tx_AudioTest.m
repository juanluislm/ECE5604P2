clear
close all
%% Parameters
% Sampling frequency
Fs = 44.1e3; 
% Duration of each sinusoid
d = 1; 
df = 0.005;
% Array of frequencies
freq_arr = (df:df:0.5).*Fs; 
Nf = length(freq_arr) % Length of frequency array
% Time domain vector
nT = (1:d*Fs)/Fs;
N_samples = length(nT); % Length of time domain vector


%% Create the sampled sines and cosines
sines = zeros(Nf, N_samples);
cosines = zeros(Nf, N_samples);
for i = 1:1:Nf
    sines(i,:) = sin(2*pi*freq_arr(i).*nT);
    cosines(i,:) = cos(2*pi*freq_arr(i).*nT);
end

for i = 1:1:Nf
    freq_arr(i)
    sound(cosines(i,:),Fs);
    pause(d)
end

