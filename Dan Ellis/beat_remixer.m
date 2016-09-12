function out = beat_remixer(filename1,filename2,subbeatlevel,mixprob,numoutbeats)

% out = beat_remixer(filename1,filename2,mixrate)
if nargin<4
    subbeatlevel = 1;
else
    if or( (subbeatlevel ~= round(subbeatlevel)), (subbeatlevel<=0) ),
        disp('please choose a whole number greater than 0 for the sub-beat-level, e.g. 2 or 4');
        out = [];
        return
    end
end

if nargin<4
    mixprob = 0.5;
end

if nargin<5
    numoutbeats =16;% number of beats for the output
end

% load the first wavefile and then extract the beats
[x1 fs] = wavread(filename1); % read in the input audio file
x1 = mean(x1,2); % convert to mono if necessary
% in case x is not at 44khz we can resample it
if (fs~=44100),
    x1 = resample(x1,44100,fs);
    fs = 44100;
end

beats1 = beat_tracker(filename1);
beats1 = interp1(1:length(beats1),beats1,1:1/subbeatlevel:length(beats1),'linear'); % sub-divide the beats according to input parameter
beatsamples1 = round(fs*beats1); % put the beat times into audio samples
beatsampledurations1 = diff(beatsamples1);
beatperiod1 = median((beatsampledurations1));

numbeatsperbar1 = 4 * subbeatlevel; % adjust for the subdivision used

% assign labels to each metrical position
beatlabels1 = 0*beats1;
for i=1:numbeatsperbar1,
    beatlabels1(i:numbeatsperbar1:end)=i;
end

lastendofbar1 = find(beatlabels1==min(beatlabels1));
lastendindex1 = lastendofbar1(end);

beats1 = beats1(1:lastendindex1);
beatsamples1 = beatsamples1(1:lastendindex1);

for k=1:length(beats1)-1,
    if mod(k,numbeatsperbar1)
        barnum1 = 1+floor(k/numbeatsperbar1);
    else
        barnum1= floor(k/numbeatsperbar1);
    end
    beatgroup1{barnum1,beatlabels1(k)}=x1(beatsamples1(k)+1:beatsamples1(k)+beatsampledurations1(k));
end


% load the second wavefile and then extract the beats
[x2 fs] = wavread(filename2); % read in the input audio file
x2 = mean(x2,2); % convert to mono if necessary
% in case x is not at 44khz we can resample it
if (fs~=44100),
    x2 = resample(x2,44100,fs);
    fs = 44100;
end

beats2 = beat_tracker(filename2);
beats2 = interp1(1:length(beats2),beats2,1:1/subbeatlevel:length(beats2),'linear'); % sub-divide the beats according to input parameter
beatsamples2 = round(fs*beats2); % put the beat times into audio samples
beatsampledurations2 = diff(beatsamples2);
beatperiod2 = median((beatsampledurations2));

numbeatsperbar2 = 4 * subbeatlevel; % adjust for the subdivision used

% assign labels to each metrical position
beatlabels2 = 0*beats2;
for i=1:numbeatsperbar2,
    beatlabels2(i:numbeatsperbar2:end)=i;
end

lastendofbar2 = find(beatlabels2==min(beatlabels2));
lastendindex2 = lastendofbar2(end);

beats2 = beats2(1:lastendindex2);
beatsamples2 = beatsamples2(1:lastendindex2);

for k=1:length(beats2)-1,
    if mod(k,numbeatsperbar2)
        barnum2 = 1+floor(k/numbeatsperbar2);
    else
        barnum2= floor(k/numbeatsperbar2);
    end
    beatgroup2{barnum2,beatlabels2(k)}=x2(beatsamples2(k)+1:beatsamples2(k)+beatsampledurations2(k));
end

% given both beat periods make a decision over whether to go faster or
% slower
% let's try slower
if (abs(1-(beatperiod1/beatperiod2))<0.05),
    stretchfact = [1 1]; % i.e. do nothing
else if beatperiod1<beatperiod2
        stretchfact = [1 (beatperiod1/beatperiod2)];
    else
        stretchfact = [(beatperiod1/beatperiod2) 1];
    end
end


out = [];
cflen = 20; % length of the cross-fade

minbeatlen = min(length(beats1),length(beats2));
% now for each beat construct a new output picking a random beat each time
for k=1:(numoutbeats*subbeatlevel),
    
    whichsong = 1+double(rand>mixprob); % trivial way to make a random choice between song 1 and 2.
    
    if (whichsong==1)
        
        % find which metrical position to choose from
        whichmetpos = mod(k,size(beatgroup1,2));
        whichmetpos(whichmetpos==0)=size(beatgroup1,2);
        
        % now determine which of the options to choose from
        whichbar = ceil(rand*size(beatgroup1,1));
        randbeat1 = beatgroup1{whichbar,whichmetpos};
        [randbeat1]=timestretch(randbeat1(:),stretchfact(1));
        
        randbeat1 = randbeat1(:)';
        
        
        % if it's the first beat, we don't need a cross-fade
        if k==1,
            out = [out randbeat1];
        else % otherwise implement a simple cross-fade over (cflen+1) samples
            out(end-cflen:end) = out(end-cflen:end).*((cflen:-1:0)/cflen);
            randbeat1(1:(cflen+1)) = randbeat1(1:(cflen+1)).*((0:1:cflen)/cflen);
            % here is the cross-fade
            out(end-cflen:end) = out(end-cflen:end) + randbeat1(1:(cflen+1));
            % now add the rest
            out = [out randbeat1((cflen+2):end)];
        end
        
    else
        % find which metrical position to choose from
        whichmetpos = mod(k,size(beatgroup2,2));
        whichmetpos(whichmetpos==0)=size(beatgroup2,2);
        
        % now determine which of the options to choose from
        whichbar = ceil(rand*size(beatgroup2,1));
        randbeat2 = beatgroup2{whichbar,whichmetpos};
        [randbeat2]=timestretch(randbeat2(:),stretchfact(2));
        randbeat2 = randbeat2(:)';
        
        % if it's the first beat, we don't need a cross-fade
        if k==1,
            out = [out randbeat2];
        else % otherwise implement a simple cross-fade over (cflen+1) samples
            out(end-cflen:end) = out(end-cflen:end).*((cflen:-1:0)/cflen);
            randbeat2(1:(cflen+1)) = randbeat2(1:(cflen+1)).*((0:1:cflen)/cflen);
            % here is the cross-fade
            out(end-cflen:end) = out(end-cflen:end) + randbeat2(1:(cflen+1));
            % now add the rest
            out = [out randbeat2((cflen+2):end)];
        end
        
        
    end
        fprintf('%.1f percent done\n',100*k/min(minbeatlen,(numoutbeats*subbeatlevel)));
    
end
    disp('100 percent done');


