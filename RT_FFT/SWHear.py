import pyaudio
import time
import numpy as np
import threading
import onset_timer
import playbeep
from scipy import interpolate
import getperiod
#from detect_peaks import detect_peaks
import getglobalcost
from scipy.signal import find_peaks_cwt
import matplotlib.pyplot as plt
class SWHear(object):

    def __init__(self,device=None,rate=None,chunk=1024):
        """fire up the SWHear class."""
        self.p=pyaudio.PyAudio()
        self.chunk = chunk #2048 # number of data points to read at a time
        self.device=device
        self.rate=rate
        self.adaptval=15
        self.conver=chunk #Length of step size...
        self.tdf=[]
        self.time_stamps=[]
        self.acorrwin=3*44100/chunk
        self.bpm=[]
        self.lentdf=300*44100/chunk
        self.thresh=0.2
        self.C=[]#np.zeros((self.acorrwin,1))
        self.fs=44100
        self.Start=[]
        self.lastbeat=[]
        self.firstrunflag = 1
        self.BT=[]
        self.acorr=[]
        self.TV=True
        self.samples=2048
        self.rayparam=np.round(43*(512/512))
        self.abpm=[]
        self.bepm = []


        self.bepmdebug=[]

    def valid_low_rate(self,device):
        """set the rate to the lowest supported audio rate."""
        for testrate in [44100]:
            if self.valid_test(device,testrate):
                return testrate
        print("SOMETHING'S WRONG! I can't figure out how to use DEV",device)
        return None

    def valid_test(self,device,rate=44100):
        """given a device ID and a rate, return TRUE/False if it's valid."""
        try:
            self.info=self.p.get_device_info_by_index(device)
            if not self.info["maxInputChannels"]>0:
                return False
            stream=self.p.open(format=pyaudio.paInt16,channels=1,
               input_device_index=device,frames_per_buffer=self.chunk,
               rate=int(self.info["defaultSampleRate"]),input=True)
            stream.close()
            return True
        except:
            return False

    def valid_input_devices(self):
        mics=[]
        for device in range(self.p.get_device_count()):
            if self.valid_test(device):
                mics.append(device)
        if len(mics)==0:
            print("no microphone devices found!")
        else:
            print("found %d microphone devices: %s"%(len(mics),mics))
        return mics

    def initiate(self):
        """run this after changing settings (like rate) before recording"""
        if self.device is None:
            self.device=self.valid_input_devices()[0] #pick the first one
        if self.rate is None:
            self.rate=self.valid_low_rate(self.device)
        if not self.valid_test(self.device,self.rate):
            print("guessing a valid microphone device/rate...")
            self.device=self.valid_input_devices()[0] #pick the first one
            self.rate=self.valid_low_rate(self.device)
        self.datax=np.arange(self.chunk)/float(self.rate)
        msg='recording from "%s" '%self.info["name"]
        msg+='(device %d) '%self.device
        msg+='at %d Hz'%self.rate
        print(msg)

    def close(self):
        """gently detach from things."""
        print(" -- sending stream termination command...")
        self.keepRecording=False #the threads should self-close
        while(self.t.isAlive()): #wait for all threads to close
            time.sleep(.1)
        self.stream.stop_stream()
        self.p.terminate()

    def stream_readchunk(self):
        """"reads some audio and re-launches itself"""
        try:
            self.data = np.fromstring(self.stream.read(self.chunk),dtype=np.int16)/1024
        except Exception as E:
            print(" -- exception! terminating...")
            print(E,"\n"*5)
            self.keepRecording=False
        if self.keepRecording:
            self.stream_thread_new()
        else:
            self.stream.close()
            self.p.terminate()
            print(" -- stream STOPPED")

    def stream_thread_new(self):
        self.t=threading.Thread(target=self.stream_readchunk)
         # self.t.daemon=True
        self.t.start()

    def stream_thread_onset(self):
        self.t2 = threading.Thread(target=self.rt_onset)
        self.t2.start()

    def stream_thread_start_beatseq(self):
        self.t4 = threading.Thread(target=self.printbeat)
        self.t4.start()

    def stream_start(self):
        """adds data to self.data until termination signal"""
        self.initiate()
        print(" -- starting stream")
        self.keepRecording=True # set this to False later to terminate stream
        self.data=None # will fill up with threaded recording data
        self.fft=None
        self.stream=self.p.open(format=pyaudio.paInt16,channels=1,
                      rate=self.rate,input=True,frames_per_buffer=self.chunk)
        self.stream_thread_new()
        self.stream_thread_onset()

    def rt_onset(self):
        self.start = time.time()
        o_step = 1024
        o_win_len = o_step * 2
        hlf_win = np.int(o_win_len / 2)
        prev_data=[]
        current_data=[]
        time.sleep(0.1)
        self.Start = time.time()
        for i in range(1,int(round(self.samples/self.chunk))):
            prev_data = np.append(prev_data,self.data)
            time.sleep(self.chunk/self.rate)

        theta1 = np.zeros(hlf_win)
        theta2 = theta1
        oldmag = theta1
        df = []
        ts=[]
        timdif=0
        count=0;
        xold=0;
        while self.TV:
            begn=time.time()
            current_data=self.data
            temp, theta1, theta2, oldmag = onset_timer.onset_detection(np.array(np.append(prev_data, current_data)),xold,
                                                                       theta1, theta2, oldmag, self.rate)
            progend=time.time()
            df = df + [temp]
            xold=current_data[len(current_data)-1]
            try:
                time.sleep(self.chunk/self.rate-progend+begn)
            except:
                count=count+1
                print(count)
                print('Error: Algo took too Long...!!!')

            prev_data = np.append(prev_data[self.chunk:len(prev_data)],current_data)
            if(len(df)>=self.adaptval):
                if(len(self.tdf)>=self.lentdf):
                    self.tdf=self.tdf[self.adaptval:len(self.tdf)]
                    self.time_stamps=self.time_stamps[self.adaptval:len(self.time_stamps)]
                df = np.array(onset_timer.part_adapt_thresh(df))
                self.time_stamps = np.concatenate((self.time_stamps, ts), axis=0)
                if (len(self.tdf) >= self.acorrwin):

                    self.pmax = np.round(60/60*(44100/self.conver))
                    self.pmin = np.round(60/120*(44100/self.conver))
                    temp = self.getbpm()
                    temp = np.array(onset_timer.adapt_threshold(list(temp)))
                    self.abpm = np.concatenate((self.abpm, temp), axis=0);
                    temp = np.argmax(temp)
                    self.bepm = self.bepm + [temp]
                    self.bpm=np.median(self.bepm[min(0,len(self.bepm)-10):])

                    self.bepmdebug = self.bepmdebug + [self.bpm]

                    self.assigncost(df)
                self.tdf = np.concatenate((self.tdf, df), axis=0)
                df=[]
                ts = []
                if timdif!=0:
                    self.conver = timdif*44100

            endn = time.time()
            ts = ts + [endn]
            timdif = endn - begn

    def assigncost(self, df):
        if self.bpm != 0 :
            if(len(self.C)>=self.lentdf):
                self.C=self.C[self.adaptval:len(self.C)]
            for j in range(0, len(df)):
                maxcost = float('-inf')
                Start = int(round(len(self.tdf) - self.bpm - self.thresh * self.bpm))
                Stop = int(round(len(self.tdf) - self.bpm + self.thresh * self.bpm+1))
                for i in range(Start, Stop):
                    if (i < 0):
                        cost = df[j]
                    else:
                        if (i > len(self.C) - 1):
                            self.C = self.tdf
                        cost = 0.9*self.C[i] + df[j]
                    if cost > maxcost:
                        maxcost = cost
                self.C = np.append(self.C, [maxcost], axis=0)
            self.lastbeat = len(self.C)-self.bpm-1+np.argmax(self.C[len(self.C)-1-self.bpm:len(self.C)-1])
            if(self.firstrunflag==1):
                self.stream_thread_start_beatseq()
                self.firstrunflag=0

    def printbeat(self):
        while self.TV:
            present_time=time.time()
            try:
                #print('BEAT')
                time.sleep(self.bpm*self.conver/self.rate-(present_time-self.time_stamps[self.lastbeat]))
                self.BT = np.append(self.BT, [time.time() - self.Start], axis=0)
                playbeep.playbeep(44100,1000,0.15)
                time.sleep(self.bpm*self.conver/self.rate/2)
            except:
                try:
                    time.sleep(2*self.bpm * self.conver / self.rate - (present_time - self.time_stamps[self.lastbeat]))
                    self.BT = np.append(self.BT, [time.time() - self.Start], axis=0)
                    playbeep.playbeep(44100, 5000, 0.15)
                    time.sleep(self.bpm * self.conver / self.rate / 2)
                except:
                    #print()
                    playbeep.playbeep(44100, 10000, 0.15)
                    time.sleep(self.bpm*self.conver/self.rate/4)

    def getbpm(self):
        n = np.arange(1,self.acorrwin)
        wv = (np.divide(n,np.power(self.rayparam,2))*np.exp(-np.divide(np.power(n,2),(2*np.power(self.rayparam ,2)))))
        eps = np.finfo(float).eps
        wv = wv / np.sum(eps + wv)

        a = self.tdf[len(self.tdf) - self.acorrwin + 1:len(self.tdf)]
        self.acorr = np.correlate(a,a, "full")
        self.acorr = self.acorr[np.int(len(self.acorr) / 2):len(self.acorr)]
        self.acorr=np.append(0,self.acorr)
        rcf = np.zeros(len(self.acorr))
        numelem = 4
        for i in range(int(self.pmin), int(self.pmax)):  # maximum beat period
            for a in range(1, numelem+1):  # number of comb elements
                for b in range(1-a,a):  # gs using normalization of comb elements
                    if((a*i+b)<len(self.acorr)):
                        rcf[i]=rcf[i]+(self.acorr[a*i+b] * wv[i])/(2*a-1)
        return rcf

if __name__=="__main__":
    ear = SWHear()
    plt.ion()
    ear.stream_start()  # goes forever
    while ear.TV:
        time.sleep(0.001)
        typedString = input()
        if typedString =='f' or typedString =='S':
            f = open('BT/DF %s'%time.time(), 'w')
            f.write("".join(str(x)+"\n" for x in ear.tdf))  # python will convert \n to os.linesep
            f2 = open('BT/tabpm %s' % time.time(), 'w')
            f2.write("".join(str(x) + "\n" for x in ear.C))  # python will convert \n to os.linesep
            f3 = open('BT/tbpm %s' % time.time(), 'w')
            f3.write("".join(str(x) + "\n" for x in ear.bepmdebug))  # python will convert \n to os.linesep
            f4 = open('BT/tbt %s' % time.time(), 'w')
            f4.write("".join(str(x) + "\n" for x in ear.BT))  # python will convert \n to os.linesep
            f2.close()
            f3.close()
            f4.close()
            f.close()

            print("Written to file!")
        if typedString =='i' or typedString =='S':
            stop=time.time()
            print(len(ear.tdf))
            ear.TV=False
            print(stop-ear.start)
            print("Interrupted By ME :P")
            ear.close()
    print("DONE")