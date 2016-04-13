clc;
close all;
load Rx_not_merged
load Rx_merged
load 6syms_test

rxs_not_merged = Rx_not_merged'.*exp(-2*pi*1i*4500/44100*[1:length(Rx_not_merged)]);
% figure (1)
% plot(rxs_not_merged)

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
newb_rxs= b_rxs_not_merged.*( -8.8079 + 4.6350i); % manually changed from "-6.0979e+02 + 4.6350e+02i" to achieve BER=0
%freqz(b_rxs_not_merged,1,2^18,44.1e3)
b_rxs = filter(shape,1,newb_rxs);
% freqz(b_rxs,1,2^18,44.1e3)
% 
% id_x = 1000:20000;
% figure(2)
% subplot(2,1,1);
% plot(real(b_rxs(id_x)));
% subplot(2,1,2);
% plot(imag(b_rxs(id_x)));

% Find Slicing point
slice_ids = (41272+(8*28)+1):k:length(b_rxs);%(m*k+(b_ord/2)+1):k:length(b_rxs) ;%(2*m*k+1):NSPS:length(frx1)
slice_ids1 = slice_ids(1:20000) ; % limit to comparison with number of symbols sent
b_rxs1 = b_rxs(slice_ids1) ;
sym1_est = sign(real(b_rxs1)) + 1i*sign(imag(b_rxs1)) ;
err_ids = find(sym1_est~=syms1) ;

err11_ids = [] ;
for nn=1:length(err_ids),
    if syms1(err_ids(nn))==1+j, err11_ids = [err11_ids err_ids(nn)] ; end ;
end ;

figure(3)
plot(real(b_rxs1),imag(b_rxs1),'x',real(b_rxs1(err11_ids)),imag(b_rxs1(err11_ids)),'ro') ;

% BER
BER = ((sum(abs(sign(real(b_rxs1(err_ids)))-real(syms1(err_ids))) ...
   + abs(sign(imag(b_rxs1(err_ids)))-imag(syms1(err_ids)))))/2)/(2*20000)