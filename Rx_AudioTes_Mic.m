clear 
close all
%%
Fs = 44100;
r = audiorecorder(Fs, 16, 1);
record(r);     % speak into microphone...

% stop(r);
%y = getaudiodata(r); 
