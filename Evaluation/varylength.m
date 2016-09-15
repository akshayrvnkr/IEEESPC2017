fs=44100;
L=6*44100:44100:30*44100;
ctr=1;
avg=[];
for l=L
    l
    runellis('/media/chaithya/Studies Related/Projects/SPC-2017/training_set/open',1,l);
    avg(ctr)=batcheval('../TestResults/GroundTruth','../TestResults/Open',l/44100);
    ctr=ctr+1;
end;