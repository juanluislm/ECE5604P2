clear 
close all
%%
Fs = 44100;
r = audiorecorder(Fs, 16, 1);
record(r);     % speak into microphone...
buff = ones(1,50000);
pause(0.1)
while(1)
    y = getaudiodata(r);
    %current = r.CurrentSample;
    buff = [buff(2:end),y(end)];

    if mean(abs(buff)) < 0.05 
        stop(r);
        %y = getaudiodata(r); 
        break;
    end
    %pause(1/Fs);
end;
