function [pos,bpm] = onset_autocorr(x,fs)
% files=dir('/media/chaithya/Studies Related/Projects/SPC-2017/training_set/open/*.wav');
% for i=1:length(files)
%     [x,fs]=audioread(sprintf('/media/chaithya/Studies Related/Projects/SPC-2017/training_set/open/%s',files(i).name));
% function to calculate the following onset detection function
% df = complex spectral difference
p=bt_parms;
% onset analysis step increment`
x=x(1:fs*2);
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
    X = fftshift(fft(win.*segment));
    
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
df=adapt_thresh(df);
%df=spline(1:length(x)/length(df):length(x),df,1:length(x));
conver=length(x)/length(df);
acorr=xcorr(df);
acorr=acorr(round(length(acorr)/2):end);
beatrange=round(44100/(3*conver)):round(44100/conver);
[pks,locs] = findpeaks(acorr(beatrange),'SortStr','descend');
bpm=44100/(3*conver)+locs(1)-1;
pos=getglobalcost(df,bpm);
pos=round(pos*conver);
end

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
end

function phase=princarg(phasein)
%phase=princarg(phasein) maps phasein into the [-pi:pi] range
phase=mod(phasein+pi,-2*pi)+pi;
end

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
end