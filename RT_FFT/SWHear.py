import pyaudio
import time
import numpy as np
import threading
import onset_timer
from detect_peaks import detect_peaks
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
        self.adaptval=15;       #TODO error for not Equal to 15 fix
        self.conver=1024; #Length of step size...
        self.tdf=[];
        self.acorrwin=85;
        self.bpm=[];
        self.fs=44100
    ### SYSTEM TESTS

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

    ### SETUP AND SHUTDOWN

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

    ### STREAM HANDLING

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
        o_step = 1024  # 1024
        o_win_len = o_step * 2
        hlf_win = np.int(o_win_len / 2)
        time.sleep(0.400)
        prev_data = self.data
        time.sleep(0.030)
        while True:
            start = time.time()
            theta1 = np.zeros(hlf_win)
            theta2 = theta1
            oldmag = theta1
            df = []
            while True:
                current_data = self.data
                temp, theta1, theta2, oldmag = onset_timer.onset_detection(np.array(np.append(prev_data, current_data)),
                                                                           theta1, theta2, oldmag, self.rate)
                df = df + [temp]
                time.sleep(0.001)
                prev_data = current_data
                if(len(df)>=self.adaptval):
                    if(len(self.tdf)>=800):
                        self.tdf=self.tdf[self.adaptval+1:len(self.tdf)];
                    df = np.array(onset_timer.adapt_threshold(df))
                    self.tdf=np.concatenate((self.tdf,df),axis=0);
                    df=[];

    def getbpm(self):
        while True:
            if(len(self.tdf)>self.acorrwin):
                acorr = np.correlate(self.tdf[len(self.tdf)-self.acorrwin+1:len(self.tdf)],self.tdf[len(self.tdf)-self.acorrwin+1:len(self.tdf)],"full")
                acorr = acorr[np.int(len(acorr)/2):len(acorr)]
                #beatrange = np.arange(self.fs/3,self.fs)/ self.conver
                #print(beatrange)
                peaks = detect_peaks(acorr(beatrange),mph=0, mpd=4)# for info look as detect_peaks
                #plt.clf()
                #plt.plot(range(np.size(acorr)),acorr)
                #plt.scatter(peaks,acorr[peaks])
                #plt.show()
                #plt.pause(0.001)
                #plt.draw()

if __name__=="__main__":
    ear = SWHear()
    plt.ion()
    ear.stream_start()  # goes forever
    while True:
        if ear.data!=None:
            plt.plot(ear.data)
            plt.show()
            plt.pause(0.1)
            plt.clf()
            plt.draw()
    print("DONE")