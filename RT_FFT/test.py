import pyaudio
import wave
import sys
import numpy as np

CHUNK = 1024

filename='/media/chaithya/Studies Related/Projects/SPC-2017/training_set/closed/closed_001.wav'
wf = wave.open(filename, 'rb')

# instantiate PyAudio (1)
p = pyaudio.PyAudio()

# open stream (2)
stream = p.open(format=p.get_format_from_width(wf.getsampwidth()),
                channels=wf.getnchannels(),
                rate=wf.getframerate(),
                input=True)

data = wf.readframes(CHUNK)

# play stream (3)
while len(data) > 0:
    print(np.fromstring(stream.read(CHUNK), dtype=np.int16) / 1024)
    #data = wf.readframes(CHUNK)


stream.stop_stream()
stream.close()

p.terminate()