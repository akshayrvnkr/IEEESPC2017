import pyaudio
import time
import numpy as np
import threading
import onset_timer
from detect_peaks import detect_peaks
import getglobalcost
from scipy.signal import find_peaks_cwt
import matplotlib.pyplot as plt
class SWHear(object):
    """
    The SWHear class is made to provide access to continuously recorded
    (and mathematically processed) microphone data.
    """
    def __init__(self,device=None,rate=None,chunk=1024):
        """fire up the SWHear class."""
        self.p=pyaudio.PyAudio()
        self.chunk = chunk #2048 # number of data points to read at a time
        self.device=device
        self.rate=rate
        self.adaptval=15       #TODO error for not Equal to 15 fix
        self.conver=1024 #Length of step size...
        self.tdf=[]
        self.time_stamps=[]
        self.acorrwin=100
        self.bpm=0
        self.lentdf=800
        self.thresh=0.2
        self.C=[]#np.zeros((self.acorrwin,1))
        self.fs=44100
        self.Start=[]
        self.lastbeat=[]
        self.firstrunflag = 1
        self.BT=[]

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
        """
        See which devices can be opened for microphone input.
        call this when no PyAudio object is loaded.
        """
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

    def stream_thread_getbpm(self):
        self.t3 = threading.Thread(target=self.getbpm)
        self.t3.start()

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
        self.stream_thread_getbpm()

    def rt_onset(self):
        start = time.time()
        o_step = 1024
        o_win_len = o_step * 2
        hlf_win = np.int(o_win_len / 2)
        time.sleep(0.400)
        prev_data = self.data
        time.sleep(self.chunk/self.rate)
        theta1 = np.zeros(hlf_win)
        theta2 = theta1
        oldmag = theta1
        df = []
        ts=[]
        self.Start = time.time()
        timdif=0
        while True:
            begn=time.time()
            current_data = self.data
            temp, theta1, theta2, oldmag = onset_timer.onset_detection(np.array(np.append(prev_data, current_data)),
                                                                       theta1, theta2, oldmag, self.rate)
            progend=time.time()
            df = df + [temp]
            try:
                time.sleep(self.chunk/self.rate-progend+begn)
            except:
                print('Error: Algo took too Long...!!!')
                #TODO what to todo???... very rare case
            prev_data = current_data
            if(len(df)>=self.adaptval):
                if(len(self.tdf)>=self.lentdf):
                    self.tdf=self.tdf[self.adaptval+1:len(self.tdf)]
                    self.time_stamps=self.time_stamps[self.adaptval+1:len(self.time_stamps)]
                df = np.array(onset_timer.adapt_threshold(df))
                self.time_stamps = np.concatenate((self.time_stamps, ts), axis=0)
                self.assigncost(df)
                self.tdf=np.concatenate((self.tdf,df),axis=0)
                df=[]
                ts=[]
                if timdif!=0:
                    self.conver = timdif*44100
            endn = time.time()
            ts = ts + [endn]
            timdif = endn - begn

    def assigncost(self, df):
        if self.bpm != 0:
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
                if (len(self.C) >= len(self.tdf)+len(df)):
                    self.C = self.C[1:len(self.C) - 1]
                self.C = np.append(self.C, [maxcost], axis=0)
            self.lastbeat = len(self.C)-self.bpm-1+np.argmax(self.C[len(self.C)-1-self.bpm:len(self.C)-1])
            if(self.firstrunflag==1):
                self.stream_thread_start_beatseq()
                self.firstrunflag=0

    def printbeat(self):
        while True:
            present_time=time.time()
            try:
                time.sleep(self.bpm*self.conver/self.rate-(present_time-self.time_stamps[self.lastbeat]))
                self.BT = np.append(self.BT, [time.time() - self.Start], axis=0)
                #print(self.BT)
            except:
                time.sleep(0.1)

    def getbpm(self):
        while True:
            if(len(self.tdf)<self.acorrwin):
                time.sleep(1)
            else:
                a = self.tdf[len(self.tdf) - self.acorrwin + 1:len(self.tdf)];
                acorr = np.correlate(a,a,"full")
                acorr = acorr[np.int(len(acorr)/2):len(acorr)]
                beatrange = np.arange(np.round(self.fs/3/self.conver),np.round(self.fs/self.conver), dtype='int32')
                peaks = detect_peaks(acorr[beatrange],mph=0, mpd=1)# for info look as detect_peaks
                #TODO Must choose a better way to select BPM .. (Rayleigh Windowing?)
                self.bpm=peaks[np.argmax(acorr[peaks])]+np.round(np.array(np.around(self.fs/3/self.conver),dtype='int32'))
                #time.sleep(0.5*self.bpm*self.conver/44100)

                #pos=getglobalcost.getglobalcost(self.tdf[len(self.tdf)-self.acorrwin+1:len(self.tdf)],peaks)
                #print(self.C)
                #plt.scatter(pos,a[pos])
                #plt.plot(self.data)
                #plt.plot(acorr);
                #plt.scatter(peaks,acorr[peaks]);
                # plt.plot(range(np.size(acorr)),acorr)
                # plt.scatter(peaks,acorr[peaks])

if __name__=="__main__":
    ear = SWHear()
    plt.ion()
    ear.stream_start()  # goes forever
    time.sleep(5)
    while True:
        if ear.tdf!=[]:
            print(time.time()-ear.time_stamps[ear.lastbeat])
            time.sleep(0.01)
    print("DONE")

    ''''def assigncost(self,df):
           for j in range(0,len(df)):
               mincost = float('inf')
               if self.bpm!=0:
                   Start=int(round(len(self.C)-self.bpm-self.thresh*self.bpm))
                   Stop=int(round(len(self.tdf)-self.bpm+self.thresh*self.bpm))
                   for i in range(Start,Stop):
                           if(i<0):
                               cost=0.9*sum(pow(np.array(self.tdf[0:len(self.tdf)-1]), 2))
                           else:
                               if (i > len(self.C) - 1):
                                   if self.C==[]:
                                       self.C=np.zeros((len(self.tdf),1))
                                   else:
                                       self.C=np.concatenate((self.C,np.zeros((len(self.tdf-len(self.C))),1)),axis=0)
                               cost = self.C[i] + 0.9*sum(pow(np.array(self.tdf[i:len(self.tdf)-1]), 2))
                   if cost < mincost:
                       mincost = cost
                   plt.clf()
                   plt.plot(self.C)
                   plt.pause(0.1)
                   plt.draw()
                   print(self.bpm)
                   if (len(self.C) >= self.lentdf):
                       self.C = self.C[1:len(self.C)-1]
                   self.C = np.append(self.C, [mincost], axis=0)'''''
