clc;
close all;
load Rx_not_merged
% load Rx_merged
load 6syms_test

Nsym = input('Number of Symbols (Press enter for Default 20000):') % number of symbols
if(isempty(Nsym))
    Nsym = 20000;
end

k = 28; % samples per symbol
m = 8; % spread
beta = 0.515; % beta
shape = rcosdesign(beta,m,k,'sqrt'); % RRC filter

fs = 44100;
nn = 1:length(Rx_not_merged);
ch = input('Channel # (Press enter for Default Ch1):'); % select channel
if(isempty(ch))
    ch = 1;
end
if(ch == 1) % Ch1 parameters
    fc = 4500;
    syms = syms1;
    syms_up = upsample(syms,k) ; 
    syms_up = [syms_up zeros(1,(2*k)+1000)]; 
    tx = filter(shape,1,syms_up) ;
end
if(ch == 2) % Ch2 parameters
    fc = 7500;
    syms = syms2;
    syms_up = upsample(syms,k) ; 
    syms_up = [syms_up zeros(1,(2*k)+1000)]; 
    tx = filter(shape,1,syms_up) ;
end
if(ch == 3) % Ch3 parameters
    fc = 10500;
    syms = syms3;
    syms_up = upsample(syms,k) ; 
    syms_up = [syms_up zeros(1,(2*k)+1000)]; 
    tx = filter(shape,1,syms_up) ;
end
if(ch == 4) % Ch4 parameters
    fc = 13500;
    syms = syms4;
    syms_up = upsample(syms,k) ; 
    syms_up = [syms_up zeros(1,(2*k)+1000)]; 
    tx = filter(shape,1,syms_up) ;
end
if(ch == 5) % Ch5 parameters
    fc = 16500;
    syms = syms5;
    syms_up = upsample(syms,k) ; 
    syms_up = [syms_up zeros(1,(2*k)+1000)]; 
    tx = filter(shape,1,syms_up) ;
end
if(ch == 6) % Ch6 parameters
    fc = 19500;
    syms = syms6;
    syms_up = upsample(syms,k) ; 
    syms_up = [syms_up zeros(1,(2*k)+1000)]; 
    tx = filter(shape,1,syms_up) ;
end

% decoding channel
rxs_not_merged = Rx_not_merged'.*exp(-2*pi*1i*fc/fs.*nn);
fac=12;
wl=800; %window length 800 seems the best length for filter
ol=wl/2; %overlap

% Lowpass Filter
b_ord = 2000;
b = firpm(b_ord,[0 0.0317 0.0363 0.5]*2,[1 1 0 0],[10 1]);
b_rxs_not_merged = filter(b,1,rxs_not_merged);

% use minimum samples (minimum of 280) after tx RRC filter to correlate signal before rx RRC filter
corr_samples = 1000;%280;
% Beginning end of correlation
a = tx(((m*k/2)+1):((m*k/2)+corr_samples));
acon = fliplr(a);
aconre = real(acon)-1i*imag(acon);
con = filter(aconre,1,b_rxs_not_merged);
[r_con,c_con] = max(abs(con));

% Finishing end of correlation
aend=tx((Nsym*k-corr_samples+1):Nsym*k);
aendcon=fliplr(aend);
aendcon=real(aendcon)-1i*imag(aendcon);
endcon=filter(aendcon,1,b_rxs_not_merged);
[r_endcon,c_endcon] = max(abs(endcon(600000:end)));
c_endcon = c_endcon + 600000-1;

% Check for rotation need or not
direc=sign(angle(endcon(c_endcon))-angle(con(c_con)));
if direc==-1
rot=0;
else
rot=-1;
end

% frequency offset & phase error correction factors
freq_factor = (rot*2*pi+(angle(endcon(c_endcon))-angle(con(c_con))))/(2*pi*(Nsym*k-corr_samples));
phase_factor = (con(c_con))';

% figure(10)
% subplot(2,1,1)
% plot(abs(con));
% subplot(2,1,2)
% plot(abs(endcon));

% corrected frequency error arising due to 2 different sound cards
for n = c_con-corr_samples+1:length(b_rxs_not_merged)
    b_rxs_new(n) = b_rxs_not_merged(n).*phase_factor*exp(-2*pi*1i*n*freq_factor);
end

% fixing phase error again
con2 = filter(aconre,1,b_rxs_new);
[r_con2,c_con2] = max(abs(con2));
new_phase_factor = (con2(c_con2))';

% Matched filter
b_rxs = filter(shape,1,b_rxs_new*new_phase_factor);

% Find Slicing point
slice_ids = (c_con-corr_samples-1+m*k/2):k:length(b_rxs);
slice_ids1 = slice_ids(1:Nsym) ; % limit to comparison with number of symbols sent
b_rxs1 = b_rxs(slice_ids1) ;
sym_est = sign(real(b_rxs1)) + 1i*sign(imag(b_rxs1)) ;
err_ids = find(sym_est~=syms(1:Nsym)) ;

% % where are the errors occurring, given that 1+j was sent?
% err11_ids = [] ;
% for nn4=1:length(err_ids),
%     if syms(err_ids(nn4))==1+1i
%         err11_ids = [err11_ids err_ids(nn4)] ; 
%     end
% end
% 
% % where are the errors occurring, given that -1+j was sent?
% err22_ids = [] ;
% for nn4=1:length(err_ids),
%     if syms(err_ids(nn4))==-1+1i
%         err22_ids = [err22_ids err_ids(nn4)] ; 
%     end
% end
% 
% % where are the errors occurring, given that -1-j was sent?
% err33_ids = [] ;
% for nn4=1:length(err_ids),
%     if syms(err_ids(nn4))==-1-1i
%         err33_ids = [err33_ids err_ids(nn4)] ; 
%     end
% end
% 
% % where are the errors occurring, given that 1-j was sent?
% err44_ids = [] ;
% for nn4=1:length(err_ids),
%     if syms(err_ids(nn4))==1-1i
%         err44_ids = [err44_ids err_ids(nn4)] ; 
%     end
% end

figure(4)
plot(real(b_rxs1),imag(b_rxs1),'x');%,real(b_rxs1(err11_ids)),imag(b_rxs1(err11_ids)),'ro',real(b_rxs1(err22_ids)),imag(b_rxs1(err22_ids)),'bo',real(b_rxs1(err33_ids)),imag(b_rxs1(err33_ids)),'yo',real(b_rxs1(err44_ids)),imag(b_rxs1(err44_ids)),'go') ;
grid on
axis square

% BER
BER = ((sum(abs(sign(real(b_rxs1(err_ids)))-real(syms(err_ids))) ...
   + abs(sign(imag(b_rxs1(err_ids)))-imag(syms(err_ids)))))/2)/(2*Nsym)