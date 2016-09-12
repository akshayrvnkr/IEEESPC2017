function runellis(folder,isopen)
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
    beats=beat_tracker([folder '\' files(j).name]);
    [x,fs]=audioread([folder '\' files(j).name]);
    dlmwrite(sprintf('../TestResults/%s/%s',type,[files(j).name(1:end-4) '.txt']),beats,'\r');
end