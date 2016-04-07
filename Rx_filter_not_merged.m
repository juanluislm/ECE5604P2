clc;
close all;
load Rx_not_merged
load Rx_merged

rxs_not_merged = Rx_not_merged'.*(exp(-2*pi*1i*4500/44100)*[1:length(Rx_not_merged)]);
%txs=real(txs);
plot(rxs_not_merged)

k = 28; % samples per symbol
m = 8; %spread
beta = 0.515; %beta

shape = rcosdesign(beta,m,k,'sqrt'); % RRC filter
fac=12;
wl=800; %window length 800 seems the best length for filter
ol=wl/2; %overlap

% Lowpass Filter
b_ord = 2000;
b = firpm(b_ord,[0 0.0317 0.0363 0.5]*2,[1 1 0 0 ], [10 1]);
b_rxs_not_merged = filter(b,1,[rxs_not_merged zeros(1,10)]);
freqz(b_rxs_not_merged,1,2^18,44.1e3)

b_rxs = filter(shape,1,b_rxs_not_merged);
freqz(b_rxs,1,2^18,44.1e3)


id_x = 1000:20000;

figure(1)
subplot(2,1,1);
plot(real(b_rxs(id_x)));
subplot(2,1,2);
plot(imag(b_rxs(id_x)));