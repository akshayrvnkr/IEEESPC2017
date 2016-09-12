function out = beat_swinger(filename,subbeatlevel,rate)

if nargin<2
    subbeatlevel = 2;
else
    if or( (subbeatlevel ~= round(subbeatlevel)), (subbeatlevel<=0) ),
        disp('please choose a whole number greater than 0 for the sub-beat-level, e.g. 2 or 4');
        out = [];
        return
    end
end

if nargin<3
    rate = 4/3;
end

if or( (rate <= 1/3), (rate>=5/3) ),
    disp('please choose a value for rate between 1/3 and 5/3, e.g. 4/3 is a good choice');
    out = [];
    return
end



% first load the wavefile and then extract the beats
[x fs] = wavread(filename); % read in the input audio file
x = mean(x,2); % convert to mono if necessary
% in case x is not at 44khz we can resample it
if (fs~=44100),
    x = resample(x,44100,fs);
    fs = 44100;
end

disp('run the beat tracker');
disp(' ');

beats = beat_tracker(filename);
beats = interp1(1:length(beats),beats,1:1/subbeatlevel:length(beats),'linear');

%beats = round(beats* 44100/1024) *1024/44100;

beatsamples = round(fs*beats); % put the beat times into audio samples

out = [];
cflen = 20;

% alternate long-short pattern
for k=1:length(beats)-1,
    slice = x(beatsamples(k)+1:beatsamples(k+1));
    if mod(k,2)
        [newslice]=timestretch(slice(:),rate);
        
    else
        [newslice]=timestretch(slice(:),2-rate);
    end
    
    newslice = newslice(:)';
    
    if k==1,
        out = [out newslice];
    else % otherwise implement a simple cross-fade over (cflen+1) samples
        out(end-cflen:end) = out(end-cflen:end).*((cflen:-1:0)/cflen);
        newslice(1:(cflen+1)) = newslice(1:(cflen+1)).*((0:1:cflen)/cflen);
        % here is the cross-fade
        out(end-cflen:end) = out(end-cflen:end) + newslice(1:(cflen+1));
        % now add the rest
        out = [out newslice((cflen+2):end)];
    end
    
    fprintf('%.1f percent done\n',100*k/length(beats));
    
end
    disp('100 percent done');

out = out(:);

