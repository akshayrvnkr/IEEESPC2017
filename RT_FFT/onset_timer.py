import time
import threading
import numpy as np
from scipy import interpolate

rpi_gpio_out = 0
lbo = 0
bri = 0
fs  = 44100
initial_waiting_time = 2 #in seconds
x = np.zeros(fs*initial_waiting_time)


def master_timer(lbo, bri):
	global rpi_gpio_out
	# TODO: 2
	rpi_gpio_out = 0
	time.sleep(bri)
	rpi_gpio_out = 1

def matlab_buffer(list,grp,overlap):
	z=(grp-len(list)%(grp-overlap))%grp
	list=list+[0]*z
	return [list[i:i+grp] for i in range(0,len(list)-overlap,grp-overlap)]

def adapt_threshold(df,pre=8,post=7):
	N=len(df)
	m=[]
	for i in range(0,min(post, N)-1):
		k = min(i + pre, N)
		m[i] = np.mean(df[0:k])

	if N > (post + pre):
		m = m + [np.mean(np.array(matlab_buffer(df, post + pre + 1, post + pre)),axis=1)]

	for i in range(N-pre,N-1):
		j = max(i - post, 1)-1
		m[i] = np.mean(df[j:len(df)-1])

	np.subtract(df,m,df)
	return (df>0)*df

def onset_detection(x,fs):
	o_step  	= 1024
	o_win_len   = o_step * 2
	hlf_win  	= o_win_len / 2
	win_hann 	= np.hanning(o_win_len)
	N 			= len(x)
	pin 		= 0
	pend		= N - o_win_len

	theta1 		= np.zeros(hlf_win)
	theta2		= np.zeros(hlf_win)
	oldmag		= np.zeros(hlf_win)

	k = 0
	df=[]
	while pin<pend:
		k += 1
		segment = x[pin : pin+o_win_len]
		x_fft	= np.fft.fft(win_hann*segment)
		x_fft	= np.fft.fftshift(x_fft)
		x_fft	= x_fft[np.floor(len(x_fft)/2):len(x_fft)]

		mag   = (np.absolute(x_fft))
		theta = np.angle(x_fft)
		dev   = ((theta-2*theta1+theta2 + np.pi) %  (-2 * np.pi)) + np.pi
		meas  = oldmag - (mag*np.exp(1j* dev))
		df[k] = np.sum(np.sqrt(np.power((np.real(meas)),2) + np.power((np.imag(meas)),2)))
		# % updatevectors
		theta2 = theta1
		theta1 = theta
		oldmag = mag
		# % move
		# to
		# next
		# frame
		pin = pin + o_step
	df = adapt_threshold(df)
	spl_tuple = interpolate.splrep(np.arange(1,len(x),len(x)/len(df)), df, s=0) # s = smoothing
	df =  interpolate.splev(np.arange(1,len(x)), spl_tuple, der=0)
	acorr = np.correlate(df,df,"full")
	acorr = acorr[np.floor(len(acorr)/2):len(acorr)]
	beatrange = np.arange(44100/3,44100)

	# df = spline(1:length(x) / length(df):length(x), df, 1:length(x));

timer_thread = threading.Thread(master_timer(lbo,bri))
timer_thread.start()
print("Abc",rpi_gpio_out)
