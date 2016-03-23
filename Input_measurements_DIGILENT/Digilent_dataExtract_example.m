clear all;      
pbr = 0.0444;
pbr_dB = 20*log10(1-pbr);
sbr_dB = -55.5;
sbr = 10^(sbr_dB/20);
pb_wl = 0.1555*2*pi; % rad/s
pb_wu = 0.2777*2*pi; % rad/s
sb_wl = 0.1111*2*pi; % rad/s
sb_wu = 0.3333*2*pi; % rad/s

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('F14example.mat');
load('unFiltered_outSineArr.mat'); unfiltered_outSineArr = outSineArr;
load('unFiltered_outCosArr.mat'); unfiltered_outCosArr = outCosArr;
d = 0.25;    % 0.25 second duration
f = 100e3;                                 T = 1/f;
nT = (1:d*f)/f;

number_sines = 27;
magArr = zeros(1,number_sines);
phaseArr = zeros(1,number_sines);

sineArr = zeros(length(nT),number_sines)     ; cosArr = zeros(length(nT),number_sines);
outSineArr = zeros(length(nT),number_sines)  ; outCosArr = zeros(length(nT),number_sines);
zArr    = zeros(number_sines,length(nT));

coeff = Rssc';
%% Create Session
s = daq.createSession('digilent')
%% Add IO channels
% The Digilent Analog Discovery board analog channels are labeled starting
% from '1':
% Analog Inputs:    '1+'/'1-' or '2+/2-'
% Analog Outputs:   'W1' / 'W2'.
s.addAnalogInputChannel('AD1',1, 'Voltage');
s.addAnalogOutputChannel('AD1',1, 'Voltage');
s.Channels

s.Rate = f;  % set rate of object to desired rate
%% Create sine/cosine files
% freqArr = [50,200,600,800,1200,1600,2000,2400,3000,3500,4000,4500,5000,...
%             6250,7500,8750,10000,11250,12500,13750,14000,14500,15000,16000,17000,18000,20000];
freqArr = [0:10:100,110:100:22.05e3];
% %% Send each sine and record output
for i = 1:1:length(freqArr)
   sineArr(:,i) = (9/10).* sin(freqArr(i)*2*pi.*nT);
   s.queueOutputData(sineArr(:,i));
   outSineArr(:,i) = s.startForeground();
end
% %% Send each cosine and record output
for i = 1:1:length(freqArr)
   cosArr(:,i)  = (9/10).* cos(freqArr(i)*2*pi.*nT);
   s.queueOutputData(cosArr(:,i));
   outCosArr(:,i) = s.startForeground();
end

% unfiltered_outSineArr = sineArr;
% unfiltered_outCosArr = cosArr;
% save('unfiltered_outSineArr','unfiltered_outSineArr_Digilent.mat');
% save('unfiltered_outCosArr','unfiltered_outCosArr_Digilent.mat');
% 
% return;

delay = 5000;
   %Create mag and phase from the last element of z (when n -> infinitiy)
%% Compare
Fs = 44100;

[MATLAB,f_Matlab] = freqz(UPRdesign,1,freqArr,Fs);
magFilt = 20.*log10(abs(MATLAB));
phaseFilt = angle(MATLAB).*(180/pi);

idx = (length(cosArr(:,1))/2):1:length(cosArr(:,1));
for i = 1:1:length(freqArr)
    
   outputArr_filtered = outCosArr(:,i)' + 1j.*outSineArr(:,i)';
   outputArr_unfiltered = unfiltered_outCosArr(:,i)' + 1j.*unfiltered_outSineArr(:,i)';
   inputArr = cosArr(:,i)'+1j*sineArr(:,i)';
   
   yn = outputArr_filtered - outputArr_unfiltered;
   zArr_yn(i,:) = yn./inputArr;
   temp = mean(zArr_yn(i,5000:length(zArr_yn(i,:))));
   magArr_yn(i) = 20.*log10(abs(temp));
   phaseArr_yn(i) = angle(temp).*(180/pi);
   
   %% Sine/cosine method
   zArr_filtered(i,:) = outputArr_filtered./inputArr; 
   zArr_unfiltered(i,:) = outputArr_unfiltered./inputArr; 
   % Find the mean value of the last values of the 'z' function
   % Then get the mag and phase
   temp = mean(zArr_filtered(i,5000:length(zArr_filtered(i,:)))); 
   temp2 = mean(zArr_unfiltered(i,5000:length(zArr_unfiltered(i,:)))); 
   
   magArr_filtered(i) = 20.*log10(abs(temp));
   phaseArr_filtered(i) = angle(temp).*(180/pi);
   % Magnitude and phase change due to the USB board (i.e. the A/D and D/A)
   magArr_unfiltered(i) = 20.*log10(abs(temp2));
   phaseArr_unfiltered(i) = angle(temp2).*(180/pi);
   
   magArr_difference(i) = magArr_filtered(i)-magArr_unfiltered(i);
   phaseArr_difference(i) = phaseArr_filtered(i)-phaseArr_unfiltered(i);
   
   %% Now do spectral analysis
   [INPUT, w] = freqz(inputArr(idx),1,2^15,Fs);
   [OUTPUT_filtered,w] = freqz(outputArr_filtered(idx),1,2^15,Fs);
   [OUTPUT_unfiltered,w] = freqz(outputArr_unfiltered(idx),1,2^15,Fs);
   
   index = find(abs(INPUT) == max(abs(INPUT)));
   magArr_filtered_spec(i) = 20.*log10(abs(OUTPUT_filtered(index)/INPUT(index)));
   phaseArr_filtered_spec(i) = angle(OUTPUT_filtered(index)/INPUT(index)).*(180/(pi));
   
   magArr_unfiltered_spec(i) = 20.*log10(abs(OUTPUT_unfiltered(index)/INPUT(index)));
   phaseArr_unfiltered_spec(i) = angle(OUTPUT_unfiltered(index)/INPUT(index)).*(180/(pi));
   
   
   magArr_difference_spec(i) = magArr_filtered_spec(i)-magArr_unfiltered_spec(i);
   phaseArr_difference_spec(i) = phaseArr_filtered_spec(i)-phaseArr_unfiltered_spec(i);
   %% Fix phase
   % For sine/cosine
   if abs(phaseFilt(i)-phaseArr_difference(i)) > 200
       if phaseArr_difference(i) > 0
           phaseArr_difference(i) = phaseArr_difference(i) - 360;
       else 
           phaseArr_difference(i) = phaseArr_difference(i) + 360;
       end
   end
   
   % For spectral analysis
   if abs(phaseFilt(i)-phaseArr_difference_spec(i)) > 200
       if phaseArr_difference_spec(i) > 0
           phaseArr_difference_spec(i) = phaseArr_difference_spec(i) - 360;
       else 
           phaseArr_difference_spec(i) = phaseArr_difference_spec(i) + 360;
       end
   end
end
%% Plot

[MATLAB,f_Matlab] = freqz(UPRdesign,1,2^15,Fs);
magFilt = 20.*log10(abs(MATLAB));
phaseFilt = angle(MATLAB).*(180/pi);


MATLAB_freqArr = freqz(UPRdesign,1,freqArr,Fs);
magFilt_freqArr = 20.*log10(abs(MATLAB_freqArr));
phaseFilt_freqArr = angle(MATLAB_freqArr).*(180/pi);

%% Response


figure(5)
p1 = plot(freqArr./Fs,phaseArr_filtered,'ro'); hold on;
p2 = plot(freqArr/Fs,phaseArr_filtered_spec, 'ks');
plot(f_Matlab./Fs,phaseFilt); 
p3 = plot(freqArr./Fs,phaseFilt_freqArr,'bx','LineWidth',1.5); title('Phase response')
hold off; 
legend([p1,p2,p3],'Sines/cosines','Windowing','MATLAB');
ylabel('Phase (degrees)'); xlabel('Fractional Frequency');

figure(6)
p1 = plot(freqArr./Fs,magArr_filtered,'ro'); hold on;
p2 = plot(freqArr/Fs,magArr_filtered_spec, 'ks');
plot(f_Matlab./Fs,magFilt);
p3 = plot(freqArr./Fs,magFilt_freqArr,'bx','LineWidth',1.5); title('Magnitude response')
hold off; 
legend([p1,p2,p3],'Sines/cosines','Windowing','MATLAB');
ylabel('Magnitude (dB)'); xlabel('Fractional Frequency');
axis([0,0.5,-120,2.5])


%% Codec

% figure(2)
% p1 = plot(freqArr./Fs,magArr_unfiltered,'bo'); hold on;
% % plot(f_Matlab./Fs,magFilt);
% % p3 = plot(freqArr./Fs,magFilt_freqArr,'bx','LineWidth',1.5); title('Magnitude response')
% hold off; 
% title('Magnitude response due to codec');
% legend([p1],'Response due to codec');
% ylabel('Magnitude (dB)'); xlabel('Fractional Frequency');
% 
% 
% figure(7)
% p1 = plot(freqArr./Fs,phaseArr_unfiltered,'ro'); hold on;
% %p3 = plot(freqArr./Fs,magArr_filtered+3.1516,'ks'); title('Magnitude response')
% plot(f_Matlab./Fs,phaseFilt);
% hold off; 
% title('Phase response due to codec');
% legend('Response due to codec','MATLAB');
% ylabel('Magnitude (dB)'); xlabel('Fractional Frequency');


%%
figure(3);
p1 = plot(freqArr./Fs,magArr_difference+3.3216,'ro'); hold on;
p2 = plot(freqArr/Fs,magArr_difference_spec+3.3216, 'ks');
plot(f_Matlab./Fs,magFilt);
p3 = plot(freqArr./Fs,magFilt_freqArr,'bx','LineWidth',1.5); title('Magnitude response with correction')
%PBspecbox(pb_wl/(2*pi),pb_wu/(2*pi), pbr_dB, -pbr_dB, [0,0,0]); hold on;
SBspecbox(0, sb_wl/(2*pi), sbr_dB, [0,0,0]); hold off;
SBspecbox(sb_wu/(2*pi), 0.5, sbr_dB, [0,0,0]);
hold off; 
legend([p1,p2,p3],'Sines/cosines','Windowing','MATLAB');
ylabel('Magnitude (dB)'); xlabel('Fractional Frequency');
axis([0,0.5,-120,2.5])

figure(4);
p1 = plot(freqArr./Fs,phaseArr_difference,'ro'); hold on;
p2 = plot(freqArr/Fs,phaseArr_difference_spec, 'ks');
plot(f_Matlab./Fs,phaseFilt); 
p3 = plot(freqArr./Fs,phaseFilt_freqArr,'bx','LineWidth',1.5); title('Phase response with correction')
hold off; 
legend([p1,p2,p3],'Sines/cosines','Windowing','MATLAB');
ylabel('Phase (degrees)'); xlabel('Fractional Frequency');

%%

%%
% save('magArr_unfiltered.mat','magArr_unfiltered');
% save('magArr_unfiltered_spec.mat','magArr_unfiltered_spec');
% 
% save('phaseArr_unfiltered.mat','phaseArr_unfiltered');
% save('phaseArr_unfiltered_spec.mat','phaseArr_unfiltered_spec');

   
   