function beats = beat_tracker(x,fs)
%
% function beats = beat_tracker(filename)
%
% this function calculates beat times in seconds from an input wavefile
p = bt_parms; % read in beat tracking parameters

df = onset_detection_function(x,p); % calculate the onset detection function
%disp('onset detection function calculated');

% get periodicity path
ppath = periodicity_path(df,p);
%disp('tempo estimation done');

% find beat locations
beats = dynamic_programming(df,p,ppath);
%disp('beat tracking completed');


% sub-functions, for simplicity all included in the same file
function p = bt_parms(res)
%  parameters for the beat tracker
if nargin<1
    res = 0.01161;
   % res=1;
end

p.fs = 44100;
p.timeres = round(p.fs * res);
p.winlen = round(512^2/p.timeres);
p.step = round(p.winlen/4);
p.bwinlen = 512; % always!
p.bstep = 128; % for the beat tracker!
% parameter for rayleigh distribution weighting
p.rayparam = round(43*(512/p.timeres));
% minimum and maximum periodicities for comb filterbank
p.pmax = round(120*(512/p.timeres));
p.pmin = round(4*(512/p.timeres));
p.lowest = round(21*(512/p.timeres)); % dixon upper limit of 247 bpm
%parameters for adaptive moving mean threshold
p.pre = round(7*(512/p.timeres));
% p.post = round(7*(512/p.timeres));
p.post = p.pre;
% factor for converting between beat period and tempo
p.fact = 60*p.fs/p.timeres;


function df = onset_detection_function(x,p)
% function to calculate the following onset detection function
% df = complex spectral difference

% onset analysis step increment`
o_step = p.winlen*2; % should be 1024
% onset analysis winlen
o_winlen = o_step*2; % should be 2048

hlfwin = o_winlen/2; % will be half fft size

% formulate hanningz window function
win = hanning(o_winlen);

% loop parameters
N = length(x);
pin = 0;
pend = N - o_winlen;

% vectors to store phase and magnitude calculations
theta1 = zeros(hlfwin,1);
theta2 = zeros(hlfwin,1);
oldmag = zeros(hlfwin,1);

% output onset detection function
df = [];

% df sample number
k = 0;
while pin<pend
    
    k=k+1;
    % calculate windowed fft frame
    segment = x(pin+1:pin+o_winlen);
    X = fft(fftshift(win.*segment));
    
    % discard first half of the spectrum
    X = X(floor(length(X)/2)+1:length(X),:);
    
    % find the magnitude and phase
    mag = (abs(X));
    theta = angle(X);
    
    % complexsd part
    dev=princarg(theta-2*theta1+theta2);
    meas=oldmag - (mag.*exp(1i.*dev));
    df(k) = sum(sqrt((real(meas)).^2+(imag(meas)).^2));
    
    % update vectors
    theta2 = theta1;
    theta1 = theta;
    oldmag = mag;
    
    % move to next frame
    pin = pin+o_step;
end

% now interpolate each detection function by a factor of 2 to get resolution of 11.6ms
df = interp1((0:length(df)-1)*p.timeres/p.fs,df,(0:0.5:length(df)-1)*p.timeres/p.fs,'cubic');
%df = interp1((0:length(df)-1)*p.timeres/p.fs,df,(0:length(df)/length(x):length(df)-1)*p.timeres/p.fs,'cubic');

function phase=princarg(phasein)
%phase=princarg(phasein) maps phasein into the [-pi:pi] range
phase=mod(phasein+pi,-2*pi)+pi;

function [ppath,obs] = periodicity_path(df,p)

% function to calculate the best periodicity path through time
% using viterbi decoding
step = p.step;
winlen = p.winlen;
n = 1:p.step;
% rayleigh weighting curve
wv = (n ./ p.rayparam .^ 2) .* exp(-n .^ 2 ./ (2*p.rayparam .^ 2));
% sum to unity
wv = sunity(wv);

pin = 0;
pend = length(df) - winlen;
% split df into overlapping frames
% find the autocorrelation function
% apply comb filtering and store output in a matrix 'obs'

ct = 0;
while(pin<pend)
    ct = ct+1;
    segment = adapt_thresh(df(pin+1:pin+winlen));
    acf(:,ct) = ftacf(segment(:));
    [rcf] = getperiod(acf(:,ct),wv,0,step,p.pmin,p.pmax);
    obs(:,ct) = sunity(adapt_thresh(rcf));
    pin = pin+step;
end

% make transition matrix
tmat = zeros(step);

% as a diagonal guassian
for i=28:108,
    tmat(:,i) = (normpdf2(n,i,8));
end; % this is a change for music therapy

tmat(1:28,:) = 0;
tmat(:,1:28) = 0;
tmat(108:128,:) = 0;
tmat(:,108:128) = 0;

% work out best path
[ppath] = viterbi_path(wv,tmat,obs+eps*max(max(obs))*rand(size(obs)));
% add on a few values at the end, to deal with final overlapping frames
ppath = [ppath ppath(end)*ones(1,4)];

function y = normpdf2(x,mu,sigma)
y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);


function [path,p] = viterbi_path(prior, transmat, obslik,filtflag)
% include acknowledgment

if nargin<4
    filtflag = 0;
end

if filtflag
    obslik = filter2(ones(5)/25,obslik);
end

scaled = 1;

T = size(obslik, 2);
prior = prior(:);
Q = length(prior);

delta = zeros(Q,T);
psi = zeros(Q,T);
path = zeros(1,T);
scale = ones(1,T);


t=1;
delta(:,t) = prior .* obslik(:,t);
if scaled
    [delta(:,t), n] = normalise(delta(:,t));
    scale(t) = 1/n;
end
psi(:,t) = 0; % arbitrary value, since there is no predecessor to t=1
for t=2:T
    for j=1:Q
        [delta(j,t), psi(j,t)] = max(delta(:,t-1) .* transmat(:,j));
        delta(j,t) = delta(j,t) * obslik(j,t);
    end
    if scaled
        [delta(:,t), n] = normalise(delta(:,t));
        scale(t) = 1/n;
    end
end
[p(T), path(T)] = max(delta(:,T));
for t=T-1:-1:1
    path(t) = psi(path(t+1),t+1);
end

if 0
    if scaled
        loglik = -sum(log(scale));
        %loglik = prob_path(prior, transmat, obslik, path);
    else
        loglik = log(p);
    end
end

function [M, z] = normalise(A, dim)
% include acknowledgment

if nargin < 2
    z = sum(A(:));
    % Set any zeros to one before dividing
    % This is vfalid, since c=0 => all i. A(i)=0 => the answer should be 0/1=0
    s = z + (z==0);
    M = A / s;
elseif dim==1 % normalize each column
    z = sum(A);
    s = z + (z==0);
    %M = A ./ (d'*ones(1,size(A,1)))';
    M = A ./ repmatC(s, size(A,1), 1);
else
    % Keith Battocchi - v. slow because of repmat
    z=sum(A,dim);
    s = z + (z==0);
    L=size(A,dim);
    d=length(size(A));
    v=ones(d,1);
    v(dim)=L;
    %c=repmat(s,v);
    c=repmat(s,v');
    M=A./c;
end

function out =sunity(in)

out = in/sum(eps+in);

function acf = ftacf(x)

x = x(:);
[M,N] = size(x);
X = fft(x,2^nextpow2(2*M-1));
acf = ifft(X.*conj(X));
acf = acf(1:M)./[M:-1:1]';


function [dfout,m] = adapt_thresh(df,pre,post)


df = df(:)';
fn = @mean;

if(nargin<2)
    pre = 8;
    post = 7;
end

% moving mean threshold

N=length(df);

for i=1:min(post,N)
    k=min(i+pre,N);
    m(i)=feval(fn,df(1:k));
end

if N>(post+pre)
    m=[m feval(fn,buffer(df,post+pre+1,post+pre,'nodelay'))];
end

for i=N+(1-pre:0)
    j=max(i-post,1);
    m(i)=feval(fn,df(j:end));
end

df = df-m;

dfout = (df>0).*df;

function [rcf] = getperiod(acf,wv,timesig,step,pmin,pmax)
rcf = zeros(1,step);

if(~timesig) % timesig unknown, must be general state
    numelem = 4;
    
    for i=pmin:pmax-1, % maximum beat period
        for a=1:numelem, % number of comb elements
            for b=1-a:a-1, % gs using normalization of comb elements
                rcf(i) = rcf(i) + (acf(a*i+b)*wv(i))/(2*a-1);
            end
        end
    end
    
else
    numelem = timesig; % timesig known must be context dependent state
    
    for i=pmin:pmax-1, % maximum beat period
        for a=1:numelem, % number of comb elements
            for b=1-a:a-1, % cds not normalizing comb elements
                rcf(i) = rcf(i) + acf(a*i+b)*wv(i);
            end
        end
    end
    
end


function [beats,localscore,cumscore] = dynamic_programming(df,p,ppath)

% function to calculate beat locations given the periodicity path
% based on Ellis "beat tracking using dynamic programming"
% but using a variable periodicity path and an optional expressive
% mode which can follow greater tempo changes

alpha = 0.95;
tightness = 5;

tempmat = (ppath'*ones(1,p.step))';
pd = round(tempmat(:));

mpd = round(median(ppath));
templt = exp(-0.5*((-mpd:mpd)/(mpd/32)).^2);
localscore = conv(templt,df);

localscore = localscore(round(length(templt)/2)+(1:length(df)));

backlink = zeros(1,length(localscore));
cumscore = zeros(1,length(localscore));

starting = 1;
for i = 1:length(localscore)
    
    prange = round(-2*pd(i)):-round(pd(i))/2;
    txwt = exp( -0.5*  (  (  ( tightness  )   *   log(prange/-pd(i))).^2)  );
    
    timerange = i + prange;
    
    zpad = max(0, min(1-timerange(1),length(prange)));
    scorecands = txwt .* [zeros(1,zpad),cumscore(timerange(zpad+1:end))];
    
    [vv,xx] = max(scorecands);
    
    cumscore(i) = alpha*vv + (1-alpha)*localscore(i);
    
    if starting == 1 && localscore(i) < 0.01*max(localscore);
        backlink(i) = -1;
    else
        backlink(i) = timerange(xx);
        starting = 0;
    end
    
end

align = getalignment2(localscore(end-512:end),ones(1,1*pd(end)),1*pd(end));

b = [];
b = length(localscore) - align;

while backlink(b(end)) > 0
    b = [b,backlink(b(end))];
end

% put the beat times into seconds
beats = sort(b)*512/44100;

function [alignment] = getalignment2(dfframe,phwv,period)

period = round(period);

dfframe = dfframe(end:-1:1);

% output of alignment comb filter
phcf = zeros(1,period);

numelem = floor((length(dfframe)-0)/period);

% fit as many as possible
for i=1:period
    for b = 1:numelem,
        phcf(i) = phcf(i) + dfframe(( b-1)*period+i) * phwv(i); %+ ...
    end
end


[val,alignment] = max(phcf);
[val2,bestguess] = max(phwv);

if alignment>=2,
    alignment = alignment-1;
end