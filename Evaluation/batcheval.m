function [avg,stddev,per] = batcheval(gtfol,outfol,End)
%%
%Usage 
%[avg,stddev,per]=batcheval('../TestResults/GroundTruth','../TestResults/Open');
gt=dir(sprintf('%s/*.txt',gtfol));
out=dir(sprintf('%s/*.txt',outfol));
per=[];
for j=1:length(out)
    per(j)=beat_evaluation(sprintf('%s/%s',gtfol,gt(j).name),sprintf('%s/%s',outfol,out(j).name),End);
end;
avg=mean(per);
stddev=sqrt(var(per));
end