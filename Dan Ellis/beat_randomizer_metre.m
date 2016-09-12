function out = beat_randomizer_metre(filename,subbeatlevel,numoutbeats)

if nargin<2
    subbeatlevel = 1;
else
    if or( (subbeatlevel ~= round(subbeatlevel)), (subbeatlevel<=0) ),
        disp('please choose a whole number greater than 0 for the sub-beat-level, e.g. 2 or 4');
        out = [];
        return
    end
end

if nargin<3
    numoutbeats =32;% number of beats for the output
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
beats = interp1(1:length(beats),beats,1:1/subbeatlevel:length(beats),'linear'); % sub-divide the beats according to input parameter
beatsamples = round(fs*beats); % put the beat times into audio samples
beatsampledurations = diff(beatsamples);

out = [];
cflen = 20; % length of the cross-fade
numbeatsperbar = 4 * subbeatlevel; % adjust for the subdivision used

% assign labels to each metrical position
beatlabels = 0*beats;
for i=1:numbeatsperbar,
    beatlabels(i:numbeatsperbar:end)=i;
end

lastendofbar = find(beatlabels==min(beatlabels));
lastendindex = lastendofbar(end);

beats = beats(1:lastendindex);
beatsamples = beatsamples(1:lastendindex);

for k=1:length(beats)-1,
    if mod(k,numbeatsperbar)
        barnum = 1+floor(k/numbeatsperbar);
    else
        barnum= floor(k/numbeatsperbar);
    end
    beatgroup{barnum,beatlabels(k)}=x(beatsamples(k)+1:beatsamples(k)+beatsampledurations(k));
end


disp('randomize the beat locations, but preserving metre');
disp(' ');

% now for each beat construct a new output picking a random beat each time
for k=1:(numoutbeats*subbeatlevel),
    
    % find which metrical position to choose from
    whichmetpos = mod(k,size(beatgroup,2));
    whichmetpos(whichmetpos==0)=size(beatgroup,2);
    
    % now determine which of the options to choose from
    whichbar = ceil(rand*size(beatgroup,1));
    randbeat = beatgroup{whichbar,whichmetpos};
    randbeat = randbeat(:)';
    
    % if it's the first beat, we don't need a cross-fade
    if k==1,
        out = [out randbeat];
    else % otherwise implement a simple cross-fade over (cflen+1) samples
        out(end-cflen:end) = out(end-cflen:end).*((cflen:-1:0)/cflen);
        randbeat(1:(cflen+1)) = randbeat(1:(cflen+1)).*((0:1:cflen)/cflen);
        % here is the cross-fade
        out(end-cflen:end) = out(end-cflen:end) + randbeat(1:(cflen+1));
        % now add the rest
        out = [out randbeat((cflen+2):end)];
    end
    
    
end

out = out(:);