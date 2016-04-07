clear 
close all
%%
Fs = 44.1e3;
r = audiorecorder(Fs, 16, 1);
record(r);     % speak into microphone...

% stop(r);
y = getaudiodata(r); 
%fd = (Fs/100:Fs/100:Fs);
%figure; plot(y);
%Rx_sines_after_sampler_long=y;
%save('Rx_sines_after_sampler_long.mat','Rx_sines_after_sampler_long');
% [Pxx, fx] = pwelch(y, hamming(1024), 512, 2^13, Fs);
% figure; plot(fx, 10*log10(Pxx))
% [H, W] = freqz(y,1,2^22,1);
% figure; plot(W, db(H))
% Rx_cosines_after_sampler = y;
% save('Rx_cosines_after_sampler.mat','Rx_cosines_after_sampler');
