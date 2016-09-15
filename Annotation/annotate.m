function annotate(filename,folder)
    %folder=/media/chaithya/Studies Related/Projects/SPC-2017/training_set/closed
    [x,fs]=audioread(sprintf('%s/%s',folder,filename));
    ctr=1;
    finish=false;
    tic;
    sound(x,fs);
    while ~finish
        k=getkey;
        if k~=113
            a(ctr)=toc;
            ctr=ctr+1;    
        else
            finish=true;
            clear sound;
        end
    end
    dlmwrite(sprintf('../TestResults/Closed/%s',[filename(1:end-4) '.txt']),a,'\r');
end