def onset(x=None, fs=None):
    # function to calculate the following onset detection function
    # df = complex spectral difference
    p = bt_parms
    # onset analysis step increment`
    x = x(mslice[1:fs * 2])
    o_step = p.winlen * 2# should be 1024
    # onset analysis winlen
    o_winlen = o_step * 2# should be 2048

    hlfwin = o_winlen / 2# will be half fft size

    # formulate hanningz window function
    win = hanning(o_winlen)

    # loop parameters
    N = length(x)
    pin = 0
    pend = N - o_winlen

    # vectors to store phase and magnitude calculations
    theta1 = zeros(hlfwin, 1)
    theta2 = zeros(hlfwin, 1)
    oldmag = zeros(hlfwin, 1)

    # output onset detection function
    df = mcat([])

    # df sample number
    k = 0
    while pin < pend:

        k = k + 1
        # calculate windowed fft frame
        segment = x(mslice[pin + 1:pin + o_winlen])
        X = fft(fftshift(win *elmul* segment))

        # discard first half of the spectrum
        X = X(mslice[floor(length(X) / 2) + 1:length(X)], mslice[:])

        # find the magnitude and phase
        mag = (abs(X))
        theta = angle(X)
        # complexsd part
        dev = princarg(theta - 2 * theta1 + theta2)
        df(k).lvalue = sum(sqrt((real(meas)) **elpow** 2 + (imag(meas)) **elpow** 2))
        # update vectors
        theta2 = theta1
        theta1 = theta
        oldmag = mag
        # move to next frame
        pin = pin + o_step
        end
        df = adapt_thresh(df)
        [pks, locs] = findpeaks(df, mstring('SortStr'), mstring('descend'))
        n = mslice[1:length(df)]
        for k in mslice[2:length(locs)]:
            impls = zeros(1, length(df))
            impls(locs(mslice[1:k])).lvalue = (pks(1) + pks(2)) / 2
            cost(k).lvalue = sum((df - impls) **elpow** 2)
            if (cost(k) > cost(k - 1) and k > 3):
                k = k - 1
                locs = locs(mslice[1:k])
                break
                end
                end
                locs = sort(locs)
                timeper = round(median(diff(locs)))
                maxp = locs(end)
                end