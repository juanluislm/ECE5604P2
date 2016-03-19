clear
close all
%% Load data
load Tx_sinusoids % = sines(Nf, N_samples), cosines(Nf, N_samples)
load Rx_sines
load Rx_cosines
%%
N_RxSines = length(Rx_sines);
N_RxCos = length(Rx_cosines);




return
Tx_SinesConcat = reshape(sines, 1,numel(sines));
[Xs,w] = freqz(Tx_SinesConcat,1,2^22,1);
[Ys,w] = freqz(Rx_sines,1,2^22,1);
figure;
plot(w,db(Xs)); hold on;
plot(w,db(Ys),'r'); legend('Tx','Rx')
