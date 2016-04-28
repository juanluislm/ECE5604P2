clc;
close all;
clear
load LP_filter
%load Tx_not_merged_sent; 
load Tx_not_merged;
load Rx_not_merged; rxs = Rx_not_merged; %rxs = Rx_not_merged(35863:length(txs)+35863-1); 
load Rx_merged
load 6syms_test; syms_concat = [syms1;syms2;syms3;syms4;syms5;syms6];
shape_length = 225; delay_shape=(shape_length-1)/2;


actual_del = 35863+4410

txs = txs(4410+delay_shape:length(txs)-4410);
%rxs = rxs(4410:end);
% % % Ptx = sum(abs(txs).^2);
% % % Prx = sum(abs(rxs).^2);
% % % P_ratio = Ptx/Prx;

%rxs = rxs.*20;%.*P_ratio;

% % % idx = 1:1:10000;
% % % figure
% % % plot(txs(idx)); hold on;
% % % plot(rxs(idx));
% % % return

%% Param
k = 28; % samples per symbol
m = 8; %spread
beta = 0.515; %beta


%% Correct

% Center frequencies
center_freq = 4500:3000:19500;
% Demodulate
% CHOOSE A CENTER FREQUENCY INDEX
center_idx = 3;
demod_fracfreq = center_freq(center_idx)/44100;

f = pi;
t = 5.612e5;
delta = f/t;

txs_not_merged = txs(1:length(txs)).*exp(-2*pi*1j*(demod_fracfreq).*(0:length(txs)-1));
%rxs_not_merged = rxs(1:t)'.*exp(-2*pi*1j*demod_fracfreq*(0:t-1) + 1j*(delta.*(0:t-1)));
rxs_not_merged = rxs(1:length(rxs))'.*exp(-2*pi*1j*(demod_fracfreq+0/44100)*(0:length(rxs)-1)).*5;


% Lowpass filter (% Length of 2000, loaded from workspace)
b_ord=2000;
b_txs_not_merged = filter(b,1,txs_not_merged);
b_rxs_not_merged = filter(b,1,rxs_not_merged);

% 
b_txs_not_merged = b_txs_not_merged(1001:end);
b_rxs_not_merged = b_rxs_not_merged(1001:end);

%% Correlate start samples
N_corr = k*5;

%C = xcorr(rxs,txs(1:N_corr));
correlator = fliplr(conj(b_txs_not_merged(1:N_corr))); 
C = filter(correlator,1,b_rxs_not_merged);

idx = find(abs(C) == max(abs(C)));
%phase0 = angle(conj(C(idx)));


delay_beg = idx-N_corr+1;

% Correct rxs
b_rxs_not_merged = b_rxs_not_merged(delay_beg+1:end).*conj(C(idx));

figure
plot(real(b_txs_not_merged)); hold on
title('txs and rxs without correction');
plot(real(b_rxs_not_merged));

figure
plot(abs(C)); title('Corr')

%% Now the end
% idx = C == max(abs(C));
% value = C(idx);
correlator2 = conj(fliplr(b_txs_not_merged(end-N_corr*5:end)));
C_end = filter(correlator2,1,b_rxs_not_merged);

idx_end = find(abs(C_end) == max(abs(C_end(fix(length(C_end)/2):end))));
max_C_end = C_end(idx_end);

phase_end = angle(conj(max_C_end));
b_rxs_not_merged = b_rxs_not_merged(1:idx_end);

delta = (phase_end)/length(b_rxs_not_merged);
b_rxs_not_merged = b_rxs_not_merged.*exp(1j*(delta.*(0:length(b_rxs_not_merged)-1)));

figure
plot(real(b_txs_not_merged)); hold on
title('txs and rxs without correction');
plot(real(b_rxs_not_merged));

% figure
% ax1=subplot(2,1,1);
% plot(real(txs_not_merged)); hold on;
% plot(real(rxs_not_merged));
% legend('Tx','Rx')
% 
% ax2=subplot(2,1,2);
% plot(imag(txs_not_merged)); hold on;
% plot(imag(rxs_not_merged));
% 
% linkaxes([ax1,ax2],'xy')

%%

%txs=real(txs);
% figure (1)
% plot(rxs_not_merged)

shape = rcosdesign(beta,m,k,'sqrt'); % RRC filter
fac=12;
wl=800; %window length 800 seems the best length for filter
ol=wl/2; %overlap

% Lowpass filter (% Length of 2000, loaded from workspace)
b_ord=2000;
% Receive filter
b_txs = filter(shape,1,b_txs_not_merged);
b_rxs = filter(shape,1,b_rxs_not_merged);

%freqz(b_rxs,1,2^18,44.1e3)

% In freq
[Tx_demod, w]= freqz(b_txs,1,2^15,44100);
[Rx_demod, w]= freqz(b_rxs,1,2^15,44100);
figure
plot(w,db(Tx_demod)); hold on;
plot(w,db(Rx_demod))
legend('Tx','Rx')

figure
plot(w,unwrap(angle(Tx_demod))); hold on;
plot(w,unwrap(angle(Rx_demod)));
legend('Tx','Rx')

% In time
%delay_LP = (length(b)-1)/2;
delay_match = (length(shape)-1);
delay = delay_match;

figure
ax1=subplot(2,1,1);
plot(real(b_txs)); hold on;
plot(real(b_rxs));
legend('Tx','Rx')

ax2=subplot(2,1,2);
plot(imag(b_txs)); hold on;
plot(imag(b_rxs));

linkaxes([ax1,ax2],'xy')
%axis([delay ,28*10 + delay ,-inf,inf])



figure;
plot(real(b_txs),imag(b_txs)); title('Tx')
figure
%plot(real(b_rxs(delay:delay+20000)),imag(b_rxs(delay:delay+20000))); title('Rx')
plot(real(b_rxs),imag(b_rxs)); title('Rx')

%id_x = 1000:20000;
% figure(2)
% subplot(2,1,1);
% plot(real(b_rxs(id_x)));
% subplot(2,1,2);
% plot(imag(b_rxs(id_x)));



% Find Slicing point
delay_LP = (length(b)-1)/2;
delay_match = (length(shape)-1)/2;
syms = syms_concat(center_idx,:);

%slice_ids = (m*k+(b_ord/2)+1):k:length(b_txs) ;%(2*m*k+1):NSPS:length(frx1)
slice_ids = (delay_match+1 :k: length(b_rxs)) ;%(2*m*k+1):NSPS:length(frx1)
slice_ids1 = slice_ids(1:15000) ; % limit to comparison with number of symbols sent
b_rxs1 = b_rxs(slice_ids1);
sym1_est = sign(real(b_rxs1)) + 1j*sign(imag(b_rxs1)) ;

syms = syms(1:15000);
err_ids = find(sym1_est~=syms);

% BER
BER = ((sum(abs(sign(real(b_rxs1(err_ids)))-real(syms(err_ids))) ...
   + abs(sign(imag(b_rxs1(err_ids)))-imag(syms(err_ids)))))/2)/(2*15000)

%%
%
% the_cancer = real(b_rxs(delay_LP+delay_match*2 :k: length(b_txs)));
% the_signal = real(syms);
% 
% the_cancer_ind = the_cancer > 0 ;
% the_cancer_indb = the_cancer < 0;
% 
% the_cancer(the_cancer_ind) = 1;
% the_cancer(the_cancer_indb)= -1;
% 
% the_signal_ind = the_signal > 0 ;
% the_signal_indb = the_signal < 0;
% 
% the_signal(the_signal_ind) = 1;
% the_signal(the_signal_indb)= -1;

