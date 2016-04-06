clc;
clear all;
close all;

Nsym = 20000 ; % number of symbols; larger number may take time but give better resolution
T = 1; %symbol time
% k, m, beta values are fixed, see variance.xlsx
k = 28; % samples per symbol
m = 8; %spread
beta = 0.515; %beta

shape = rcosdesign(beta,m,k,'sqrt'); % RRC filter

NSPS = k ; % samples per symbol

syms1 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; % Dirac deltas
syms2 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 
syms3 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 
syms4 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 
syms5 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 
syms6 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 

syms_up1 = upsample(syms1,NSPS) ; %upsampling
syms_up2 = upsample(syms2,NSPS) ;
syms_up3 = upsample(syms3,NSPS) ;
syms_up4 = upsample(syms4,NSPS) ;
syms_up5 = upsample(syms5,NSPS) ;
syms_up6 = upsample(syms6,NSPS) ;

syms_up1 = [syms_up1 zeros(1,(2*NSPS)+1000)]; %Add 0s
syms_up2 = [syms_up2 zeros(1,(2*NSPS)+1000)];
syms_up3 = [syms_up3 zeros(1,(2*NSPS)+1000)];
syms_up4 = [syms_up4 zeros(1,(2*NSPS)+1000)];
syms_up5 = [syms_up5 zeros(1,(2*NSPS)+1000)];
syms_up6 = [syms_up6 zeros(1,(2*NSPS)+1000)];

tx1 = filter(shape,1,syms_up1) ; % tx after filter
tx2 = filter(shape,1,syms_up2) ;
tx3 = filter(shape,1,syms_up3) ;
tx4 = filter(shape,1,syms_up4) ;
tx5 = filter(shape,1,syms_up5) ;
tx6 = filter(shape,1,syms_up6) ;

nn1 = 1:length(tx1); %length of tx
nn2 = 1:length(tx2);
nn3 = 1:length(tx3);
nn4 = 1:length(tx4);
nn5 = 1:length(tx5);
nn6 = 1:length(tx6);

tx_fc1 = tx1.*(exp(2*pi*1i*4500/44100.*nn1)); 
tx_fc2 = tx2.*(exp(2*pi*1i*7500/44100.*nn2));
tx_fc3 = tx3.*(exp(2*pi*1i*10500/44100.*nn3));
tx_fc4 = tx4.*(exp(2*pi*1i*13500/44100.*nn4));
tx_fc5 = tx5.*(exp(2*pi*1i*16500/44100.*nn5));
tx_fc6 = tx6.*(exp(2*pi*1i*19500/44100.*nn6));

tx1 = tx_fc1; %tx signal with real and imag componenet
tx2 = tx_fc2;
tx3 = tx_fc3;
tx4 = tx_fc4;
tx5 = tx_fc5;
tx6 = tx_fc6;

tx1r = real(tx_fc1); %tx signal only real component
tx2r = real(tx_fc2);
tx3r = real(tx_fc3);
tx4r = real(tx_fc4);
tx5r = real(tx_fc5);
tx6r = real(tx_fc6);

txs=tx1r+tx2r+tx3r+tx4r+tx5r+tx6r; % sum of txs

fac=12;
wl=800; %window length 800 seems the best length for filter
ol=wl/2; %overlap

[fft1,f1]=pwelch(tx1,hamming(wl),ol,2^fac,44100,'centered');
[fft2,f2]=pwelch(tx2,hamming(wl),ol,2^fac,44100,'centered');
[fft3,f3]=pwelch(tx3,hamming(wl),ol,2^fac,44100,'centered');
[fft4,f4]=pwelch(tx4,hamming(wl),ol,2^fac,44100,'centered');
[fft5,f5]=pwelch(tx5,hamming(wl),ol,2^fac,44100,'centered');
[fft6,f6]=pwelch(tx6,hamming(wl),ol,2^fac,44100,'centered');

[fft1r,f1r]=pwelch(tx1r,hamming(wl),ol,2^fac,44100,'centered');
[fft2r,f2r]=pwelch(tx2r,hamming(wl),ol,2^fac,44100,'centered');
[fft3r,f3r]=pwelch(tx3r,hamming(wl),ol,2^fac,44100,'centered');
[fft4r,f4r]=pwelch(tx4r,hamming(wl),ol,2^fac,44100,'centered');
[fft5r,f5r]=pwelch(tx5r,hamming(wl),ol,2^fac,44100,'centered');
[fft6r,f6r]=pwelch(tx6r,hamming(wl),ol,2^fac,44100,'centered');

[ffts,fts]=pwelch(txs,hamming(wl),ol,2^fac,44100,'centered'); %sum

% tx# are tx with real+imag component
% tx#r are tx with only real componenet
% txs is sum of tx#r

% To see real+imag componenet txs, plot first 6 plots
% To see real component txs, plot 7~12 plots
% To see sum of txs, plot last one

figure(1)
% plot(fftshift(f1)/1000,fftshift(db(fft1)),'b'); %PSD of tx (real+imag)
% hold on
% plot(fftshift(f2)/1000,fftshift(db(fft2)),'r'); 
% plot(fftshift(f3)/1000,fftshift(db(fft3)),'k'); 
% plot(fftshift(f4)/1000,fftshift(db(fft4)),'m'); 
% plot(fftshift(f5)/1000,fftshift(db(fft5)),'g'); 
% plot(fftshift(f6)/1000,fftshift(db(fft6)),'c'); 

plot(fftshift(f1r)/1000,fftshift(db(fft1r)),'b+'); %PSD of tx(real)
hold on
plot(fftshift(f2r)/1000,fftshift(db(fft2r)),'r+'); 
plot(fftshift(f3r)/1000,fftshift(db(fft3r)),'k+'); 
plot(fftshift(f4r)/1000,fftshift(db(fft4r)),'m+'); 
plot(fftshift(f5r)/1000,fftshift(db(fft5r)),'g+'); 
plot(fftshift(f6r)/1000,fftshift(db(fft6r)),'c+'); 

plot(fftshift(fts)/1000,fftshift(db(ffts)),'r-S'); %PSD of sum tx

grid on
xlim([0 23]) %change window xaxis
legend('N1','N2','N3','N4','N5','N6','sum')
title('beta= 0.515, m=8, k=28, Nsym=1000, T=1')
ylabel('Mag(dB)')
xlabel('f(kHz)')
