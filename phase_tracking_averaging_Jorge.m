clc;
close all;
clear
load LP_filter
my_BERs = zeros(1,20);
% CHOOSE A CENTER FREQUENCY INDEX
center_idx = 2;

%load Tx_not_merged_sent; 
load Tx_not_merged;
load Rx_not_merged; rxs = Rx_not_merged; %rxs = Rx_not_merged(35863:length(txs)+35863-1); 
load Rx_merged
load 6syms_test; syms_concat = [syms1;syms2;syms3;syms4;syms5;syms6];
load tx_MATLAB; txs_concat = [tx1;tx2;tx3;tx4;tx5;tx6];
shape_length = 225; delay_shape=(shape_length-1)/2;


actual_del = 35863+4410;

txs = txs_concat(center_idx,:);
txs = txs(delay_shape:length(txs));

%% Param
k = 28; % samples per symbol
m = 8; %spread
beta = 0.515; %beta


%% Correct

% Center frequencies
center_freq = 4500:3000:19500;
% Demodulate
demod_fracfreq = center_freq(center_idx)/44100;

f = pi;
t = 5.612e5;
delta = f/t;

rxs_not_merged = rxs(1:length(rxs))'.*exp(-2*pi*1j*(demod_fracfreq+0/44100)*(0:length(rxs)-1)).*5;


% Lowpass filter (% Length of 2000, loaded from workspace)
b_rxs_not_merged = filter(b,1,rxs_not_merged);
b_rxs_not_merged = b_rxs_not_merged(1001:end);

%% Correlate start samples
b_txs_not_merged = txs;
N_corr = k*5;
correlator = fliplr(conj(b_txs_not_merged(1:N_corr))); 
C = filter(correlator,1,b_rxs_not_merged);
idx = find(abs(C) == max(abs(C)));
%phase0 = angle(conj(C(idx)));
delay_beg = idx-N_corr+1;
% Correct rxs
b_rxs_not_merged = b_rxs_not_merged(delay_beg+1:end).*conj(C(idx));

shape = rcosdesign(beta,m,k,'sqrt'); % RRC filter
fac=12;
wl=800; %window length 800 seems the best length for filter
ol=wl/2; %overlap

% Lowpass filter (% Length of 2000, loaded from workspace)
b_ord=2000;
% Receive filter
b_txs = filter(shape,1,b_txs_not_merged);
b_rxs = filter(shape,1,b_rxs_not_merged);
b_txs = b_txs(72:end);
b_rxs = b_rxs(72:end);
% Delay of filter
delay_match = (length(shape)-1);
delay = delay_match;

% Correct for delay
%delay_LP = (length(b)-1)/2;
delay_match = (length(shape)-1);
delay = delay_match;
syms = syms_concat(center_idx,:); % ACTUAL SYMBOLS
% Find slice ID's
slice_ids = (1 :k: length(b_rxs)) ;
slice_ids1 = slice_ids(1:15000) ; % limit to comparison with number of symbols sent
b_rxs1 = b_rxs(slice_ids1); 
sym1_est = sign(real(b_rxs1)) + 1j*sign(imag(b_rxs1)) ; % <---SYMBOL ESTIMATES

syms = sign(real(syms)) + 1j*sign(imag(syms)) ;
syms = syms(1:15000); % Cut Tx signal down (for now)

counter = 0;
% for i=1:1:length(syms)
%    if (real(sym1_est(i))~=real(syms(i)))&&(imag(sym1_est(i))~=imag(syms(i))), 
%         counter = counter+1;
%         err_ids(counter) = i;
%    end;
%        
% end
err_ids = find(sym1_est~=syms);


% figure
% ax1=subplot(2,1,1);
% plot(real(b_txs)); hold on;
% plot(real(b_rxs));
% legend('Tx','Rx')
% ax2=subplot(2,1,2);
% plot(imag(b_txs)); hold on;
% plot(imag(b_rxs));
% linkaxes([ax1,ax2],'xy')
% %axis([delay ,28*10 + delay ,-inf,inf])

figure
ax1=subplot(2,1,1);
plot(real(sym1_est)); hold on;
plot(real(syms));
legend('Tx','Rx')
ax2=subplot(2,1,2);
plot(imag(sym1_est)); hold on;
plot(imag(syms));
linkaxes([ax1,ax2],'xy')
axis([-inf,inf,-1.2,1.2])
%axis([delay ,28*10 + delay ,-inf,inf]

return
%% Track the changes
N_corr_old = N_corr;
decision_points = [angle(1+1j),angle(-1+1j),angle(-1,-1j),angle(1,-1j)];
counter = 0;
for i = 1:length(b_rxs) 
    counter = counter+1;
    
    f= [1990 1998 2001 2004]
    val = 2000 %value to find
    tmp = abs(f-val)
    [idx idx] = min(tmp) %index of closest value
    closest = f(idx) %closest value
    
    
end;





return














% save('new_b_rx_252pm.mat','b_rxs_not_merged');
% figure
% plot(real(b_txs_not_merged)); hold on
% title('txs and rxs with correction');
% plot(real(b_rxs_not_merged));

% figure
% plot(abs(C)); title('Corr')

%%

%txs=real(txs);
% figure (1)
% plot(rxs_not_merged)
return

fac=12;
wl=800; %window length 800 seems the best length for filter
ol=wl/2; %overlap

% Lowpass filter (% Length of 2000, loaded from workspace)
b_ord=2000;
% Receive filter
b_txs = filter(shape,1,b_txs_not_merged);
b_rxs = filter(shape,1,b_rxs_not_merged);



return
% figure
% ax1=subplot(2,1,1);
% plot(real(b_txs)); hold on;
% plot(real(b_rxs));
% legend('Tx','Rx')
% 
% ax2=subplot(2,1,2);
% plot(imag(b_txs)); hold on;
% plot(imag(b_rxs));
% 
% linkaxes([ax1,ax2],'xy')
% %axis([delay ,28*10 + delay ,-inf,inf])
% 
% 
% 
% figure;
% plot(real(b_txs),imag(b_txs)); title('Tx')
% figure
% %plot(real(b_rxs(delay:delay+20000)),imag(b_rxs(delay:delay+20000))); title('Rx')
% plot(real(b_rxs),imag(b_rxs)); title('Rx')

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

my_BERs(tt) = BER;
