clc;
clear
close all;
load Tx_not_merged
load 6syms_test

Nsym = 20000 ; % number of symbols; larger number may take time but give better resolution
T = 1; %symbol time
% k, m, beta values are fixed, see variance.xlsx
k = 28; % samples per symbol
m = 8; %spread
beta = 0.515; %beta

shape = rcosdesign(beta,m,k,'sqrt'); % RRC filter

NSPS = k ; % samples per symbol

% syms1 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; % Dirac deltas
% syms2 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 
% syms3 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 
% syms4 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 
% syms5 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 
% syms6 = sign(randn(1,Nsym)) + 1i*sign(randn(1,Nsym)) ; 

syms_up1 = upsample(syms1,NSPS) ; %upsampling
syms_up2 = upsample(syms2,NSPS) ;
syms_up3 = upsample(syms3,NSPS) ;
syms_up4 = upsample(syms4,NSPS) ;
syms_up5 = upsample(syms5,NSPS) ;
syms_up6 = upsample(syms6,NSPS) ;

% syms_up1 = [syms_up1 zeros(1,(2*NSPS)+1000)]; %Add 0s
% syms_up2 = [syms_up2 zeros(1,(2*NSPS)+1000)];
% syms_up3 = [syms_up3 zeros(1,(2*NSPS)+1000)];
% syms_up4 = [syms_up4 zeros(1,(2*NSPS)+1000)];
% syms_up5 = [syms_up5 zeros(1,(2*NSPS)+1000)];
% syms_up6 = [syms_up6 zeros(1,(2*NSPS)+1000)];

tx1 = filter(shape,1,syms_up1) ; % tx after filter
tx2 = filter(shape,1,syms_up2) ;
tx3 = filter(shape,1,syms_up3) ;
tx4 = filter(shape,1,syms_up4) ;
tx5 = filter(shape,1,syms_up5) ;
tx6 = filter(shape,1,syms_up6) ;


save('tx_MATLAB.mat','tx1','tx2','tx3','tx4','tx5','tx6');
return

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
txs = [zeros(1,0.1*44.1e3), txs, zeros(1,0.1*44.1e3)];

%sound(txs,44.1e3)


