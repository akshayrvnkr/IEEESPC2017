function [beat_samples,odf] = start
    x = []; 
    framelen = 512;
    winlen = 512;
    odf = [];
    theta1 = zeros(winlen/2,1);
    theta2 = zeros(winlen/2,1);
    oldmag = zeros(winlen/2,1);
    beat_samples = [];
    count = 0;
    ip = audioread('120bpm.wav');
    while(length(x)<length(ip))
        count = count + 1;
        if ((length(ip)-length(x))>framelen)
            x = vertcat(x,ip((count-1)*framelen + 1:(count)*framelen));
            [odf,theta1,theta2,oldmag] = get_onset(x,odf,theta1,theta2,oldmag);
            beat_samples = analyze_frame(odf,beat_samples,44100);
            pause(0.01);
        else
            x = vertcat(x,ip((count-1)*framelen + 1:end));
            [odf,theta1,theta2,oldmag] = get_onset(x,odf,theta1,theta2,oldmag);
            break;
        end
    end
    beat_locations = zeros(1,length(odf));
    beat_locations(beat_samples) = 30;
    stem(beat_locations)
    hold on
    plot(odf)
end

function [odf,theta1,theta2,oldmag] = get_onset(x,odf,theta1,theta2,oldmag)
    framelen = 512;
    winlen = 512;
    
    % From the entire data obtained till now, choose 50 % from the current
    % frame and 50 % from previous frame
    if length(x) > framelen
        x_chunk = x(end - 2*winlen + 1:end - winlen);
    
        % Pass the chunk to obtain the onset value and update the onset
        % detection function
        
        [temp,theta1,theta2,oldmag] = onset(x_chunk,theta1,theta2,oldmag);
        odf = horzcat(odf,temp);
    end
end

function beat_samples = analyze_frame(odf,beat_samples,fs)

    Bf = 512;
    Bh = Bf / 4;
    if ((length(odf) >= Bf) && (mod(length(odf),Bh) == 0))
        ana_frame = odf(end-Bf+1:end);
        ana_frame_mean = adapt_mean_threshold(ana_frame);
        ana_frame = ((ana_frame - ana_frame_mean) + abs(ana_frame - ana_frame_mean)) / 2;
        acorr = (ifft(fft(ana_frame) .* conj(fft(ana_frame))));
        % plot(acorr)
        yg = zeros(Bh,1);
        for i = 1:Bh
            for j = 1:Bf
                yg(i) = yg(i) + (acorr(j)*fglt(j,i));
            end
        end
        % plot(yg)
        [~,taug] = max(yg);
        bpm = 60 / taug / (512 / fs);
        display(bpm)
        ana_frame = fliplr(ana_frame);
        Hg = gethg(taug,Bf);
        zg = Hg' * ana_frame';
        [~,alphag] = max(zg);
        i = 1;
        while(i > 0)
            if (alphag + (i-1)*taug) < Bh
                temp1(i) = length(odf) + alphag + (i)*taug;
                i = i + 1;
            else
                break;
            end
        end
        beat_samples = horzcat(beat_samples,temp1);
    end
end

function [df,theta1,theta2,oldmag] = onset(x_chunk,theta1,theta2,oldmag)
    framelen = 512;
    winlen = 512;
    win = hanning(winlen);
    X = fft(fftshift(win.*x_chunk));
    X = X(floor(length(X)/2)+1:length(X),:);
    mag = (abs(X));
    theta = angle(X);
    dev = princarg(theta-2*theta1+theta2);
    meas = oldmag - (mag.*exp(1i.*dev));
    df = sum(sqrt((real(meas)).^2+(imag(meas)).^2));
    theta2 = theta1;
    theta1 = theta;
    oldmag = mag;
end

function phase = princarg(phasein)
    phase = mod(phasein+pi,-2*pi)+pi;
end

function ana_frame_mean = adapt_mean_threshold(ana_frame)
    ana_frame_mean = ana_frame;
    q = 16;
    for i = 1:length(ana_frame)
        if i <= q/2
            ana_frame_mean(i) = mean(ana_frame(i:i+(q/2)));
        elseif (i > q/2) && (i < length(ana_frame)-q/2)
            ana_frame_mean(i) = mean(ana_frame(i-(q/2):i+(q/2)));
        elseif i >= length(ana_frame)-q/2
            ana_frame_mean(i) = mean(ana_frame(i-(q/2):i));
        end
    end
end

function fg = fglt(l,tau)
    beta = 43;
    lambda = 0;
    for p = 1:4
        for v = 1-p:p-1
            if l == (tau*p - v)
                lambda = 1 / (2*p - 1);
            end
        end
    end
    fg = (lambda * tau / beta / beta) * exp(- (tau*tau)/(2*beta*beta));
end

function Hg = gethg(taug,Bf)
    Hg = zeros(Bf,taug);
    for i = 1:taug
        for j = 1:Bf
            for k = 1:floor(Bf/taug)
                if j == ((k*taug) - i)
                    Hg(j,i) = (Bf - j) / Bf; 
                end
            end
        end
    end  
end