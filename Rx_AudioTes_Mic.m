clear 
close all
%%
Fs = 44.1e3;
r = audiorecorder(Fs, 16, 1);
record(r);     % speak into microphone...

% stop(r);
%y = getaudiodata(r); 
