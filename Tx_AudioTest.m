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


%% Create the sampled sines and cosines
sines = zeros(Nf, N_samples);
%cosines = zeros(Nf, Nf*Fs*(d + dd));
for i = 1:1:Nf
    sines(i,:) = sin(2*pi*freq_arr(i).*nT);
    %cosines( ((i-1)*Fs*d +(i-1)*Fs*dd +1):i*Fs*(d + dd) )  = [cos(2*pi*freq_arr(i).*nT), zeros(dd*Fs,1)];
end
[n,m] = size(sines);
sines(:,m-(dd*Fs)+1:m) = zeros(Nf,dd*Fs);
sines_col = sines';


sines = reshape(sines_col,1,numel(sines_col));

%% Send
sound(sines,Fs);
%pause(d*2)

% clear sound
