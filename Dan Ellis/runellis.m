function runellis(folder,isopen,len)
%This saves the results in TestResults/Open or Closed folder based on if
%isopen is 1 or 0
type='Open';
if(~isopen)
    type='Closed';
end;
files=dir(sprintf('%s/*.wav',folder));
for j=1:length(files)
    j
    warning('off');
    [x,fs]=audioread([folder '/' files(j).name]);
    if(len==[])
        len=length(x);
    end;
    if (fs~=44100),
        x = resample(x,44100,fs);
    end
    x=x(1:len);
    x = mean(x,2);
    beats=beat_tracker(x,fs);
    dlmwrite(sprintf('../TestResults/%s/%s',type,[files(j).name(1:end-4) '.txt']),beats,'\r');
end;
end