clearvars();
outputdir = ('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out\auditory');

varNames = {'TTP1', 'TTP2', 'TTP3', 'AMP1', 'AMP2', 'AMP3', 'ONS1', 'ONS2', 'ONS3'};
varTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
stattable{1} = table('Size', [10000, 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{2} = table('Size', [10000, 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{3} = table('Size', [10000, 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);

stattable{4} = table('Size', [10000, 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{5} = table('Size', [10000, 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{6} = table('Size', [10000, 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);

load("E:\boot_out\auditory\OutputMatrix");

for sampleidx=1:size(OutputMatrix{1,1},2)
stattable{1}(sampleidx,1:3) = table(OutputMatrix{1,1}(sampleidx).MaxIndexInSeconds, OutputMatrix{2,1}(sampleidx).MaxIndexInSeconds, OutputMatrix{3,1}(sampleidx).MaxIndexInSeconds);
stattable{1}(sampleidx,4:6) = table(OutputMatrix{1,1}(sampleidx).MaxValue, OutputMatrix{2,1}(sampleidx).MaxValue, OutputMatrix{3,1}(sampleidx).MaxValue);
stattable{1}(sampleidx, 7:9) = table(OutputMatrix{1,1}(sampleidx).OnsetTimeInSeconds, OutputMatrix{2,1}(sampleidx).OnsetTimeInSeconds, OutputMatrix{3,1}(sampleidx).OnsetTimeInSeconds);

stattable{2}(sampleidx,1:3) = table(OutputMatrix{1,2}(sampleidx).MaxIndexInSeconds, OutputMatrix{2,2}(sampleidx).MaxIndexInSeconds, OutputMatrix{3,2}(sampleidx).MaxIndexInSeconds);
stattable{2}(sampleidx,4:6) = table(OutputMatrix{1,2}(sampleidx).MaxValue, OutputMatrix{2,2}(sampleidx).MaxValue, OutputMatrix{3,2}(sampleidx).MaxValue);
stattable{2}(sampleidx, 7:9) = table(OutputMatrix{1,2}(sampleidx).OnsetTimeInSeconds, OutputMatrix{2,2}(sampleidx).OnsetTimeInSeconds, OutputMatrix{3,2}(sampleidx).OnsetTimeInSeconds);

stattable{3}(sampleidx,1:3) = table(OutputMatrix{1,3}(sampleidx).MaxIndexInSeconds, OutputMatrix{2,3}(sampleidx).MaxIndexInSeconds, OutputMatrix{3,3}(sampleidx).MaxIndexInSeconds);
stattable{3}(sampleidx,4:6) = table(OutputMatrix{1,3}(sampleidx).MaxValue, OutputMatrix{2,3}(sampleidx).MaxValue, OutputMatrix{3,3}(sampleidx).MaxValue);
stattable{3}(sampleidx, 7:9) = table(OutputMatrix{1,3}(sampleidx).OnsetTimeInSeconds, OutputMatrix{2,3}(sampleidx).OnsetTimeInSeconds, OutputMatrix{3,3}(sampleidx).OnsetTimeInSeconds);

stattable{4}(sampleidx,1:3) = table(OutputMatrix{1,1}(sampleidx).MaxIndexInSeconds, OutputMatrix{1,2}(sampleidx).MaxIndexInSeconds, OutputMatrix{1,3}(sampleidx).MaxIndexInSeconds);
stattable{4}(sampleidx,4:6) = table(OutputMatrix{1,1}(sampleidx).MaxValue, OutputMatrix{1,2}(sampleidx).MaxValue, OutputMatrix{1,3}(sampleidx).MaxValue);
stattable{4}(sampleidx, 7:9) = table(OutputMatrix{1,1}(sampleidx).OnsetTimeInSeconds, OutputMatrix{1,2}(sampleidx).OnsetTimeInSeconds, OutputMatrix{1,3}(sampleidx).OnsetTimeInSeconds);

stattable{5}(sampleidx,1:3) = table(OutputMatrix{2,1}(sampleidx).MaxIndexInSeconds, OutputMatrix{2,2}(sampleidx).MaxIndexInSeconds, OutputMatrix{2,3}(sampleidx).MaxIndexInSeconds);
stattable{5}(sampleidx,4:6) = table(OutputMatrix{2,1}(sampleidx).MaxValue, OutputMatrix{2,2}(sampleidx).MaxValue, OutputMatrix{2,3}(sampleidx).MaxValue);
stattable{5}(sampleidx, 7:9) = table(OutputMatrix{2,1}(sampleidx).OnsetTimeInSeconds, OutputMatrix{2,2}(sampleidx).OnsetTimeInSeconds, OutputMatrix{2,3}(sampleidx).OnsetTimeInSeconds);

stattable{6}(sampleidx,1:3) = table(OutputMatrix{3,1}(sampleidx).MaxIndexInSeconds, OutputMatrix{3,2}(sampleidx).MaxIndexInSeconds, OutputMatrix{3,3}(sampleidx).MaxIndexInSeconds);
stattable{6}(sampleidx,4:6) = table(OutputMatrix{3,1}(sampleidx).MaxValue, OutputMatrix{3,2}(sampleidx).MaxValue, OutputMatrix{3,3}(sampleidx).MaxValue);
stattable{6}(sampleidx, 7:9) = table(OutputMatrix{3,1}(sampleidx).OnsetTimeInSeconds, OutputMatrix{3,2}(sampleidx).OnsetTimeInSeconds, OutputMatrix{3,3}(sampleidx).OnsetTimeInSeconds);
end

TTPstattableL = [table2array(stattable{1}(:,1:3)), table2array(stattable{2}(:,1:3)), table2array(stattable{3}(:,1:3))];
TTPstattableI = [table2array(stattable{4}(:,1:3)), table2array(stattable{5}(:,1:3)), table2array(stattable{6}(:,1:3))];
AMPstattableL = [table2array(stattable{1}(:,4:6)), table2array(stattable{2}(:,4:6)), table2array(stattable{3}(:,4:6))];
AMPstattableI = [table2array(stattable{4}(:,4:6)), table2array(stattable{5}(:,4:6)), table2array(stattable{6}(:,4:6))];
ONSstattableL = [table2array(stattable{1}(:,7:9)), table2array(stattable{2}(:,7:9)), table2array(stattable{3}(:,7:9))];
ONSstattableI = [table2array(stattable{4}(:,7:9)), table2array(stattable{5}(:,7:9)), table2array(stattable{6}(:,7:9))];


comparisons = ["durations", "durations", "durations", "intensities", "intensities", "intensities"];
levels = ["intensity1", "intensity2", "intensity3", "duration1", "duration2", "duration3"];


TTP = TTPstattableL; %table2array(stattable{x}(:,1:3));
AMP = AMPstattableL; %table2array(stattable{x}(:,4:6));
ONS = ONSstattableL; %table2array(stattable{x}(:,7:9));

subjmTTP = nanmean(TTP,2);
grpmTTP = nanmean(TTP,'all');
subjmAMP = nanmean(AMP,2);
grpmAMP = nanmean(AMP,'all');
subjmONS = nanmean(ONS,2);
grpmONS = nanmean(ONS,'all'); 

for z = 1:size(TTP,1)
    TTPy(z,:) = TTP(z,:) - subjmTTP(z,:) + grpmTTP;
    AMPy(z,:) = AMP(z,:) - subjmAMP(z,:) + grpmAMP;
    ONSy(z,:) = ONS(z,:) - subjmONS(z,:) + grpmONS;
end

for y = 1:size(TTP,2)
    TTPstat(1,y) = nanmean(TTP(:,y),1); %means
    AMPstat(1,y) = nanmean(AMP(:,y),1);
    ONSstat(1,y) = nanmean(ONS(:,y),1);
    
    TTPstat(2,y) = nanstd(TTP(:,y),1); %sd
    AMPstat(2,y) = nanstd(AMP(:,y),1);
    ONSstat(2,y) = nanstd(ONS(:,y),1);
    
    TTPstat(3,y) = prctile(TTP(:,y), 2.5, 1); %ci
    TTPstat(4,y) = prctile(TTP(:,y), 97.5, 1);
    AMPstat(3,y) = prctile(AMP(:,y), 2.5, 1);
    AMPstat(4,y) = prctile(AMP(:,y), 97.5, 1);
    ONSstat(3,y) = prctile(ONS(:,y), 2.5, 1);
    ONSstat(4,y) = prctile(ONS(:,y), 97.5, 1);

    TTPstat(5,y) = prctile(TTPy(:,y), 2.5, 1); %ci-within
    TTPstat(6,y) = prctile(TTPy(:,y), 97.5, 1);
    AMPstat(5,y) = prctile(AMPy(:,y), 2.5, 1);
    AMPstat(6,y) = prctile(AMPy(:,y), 97.5, 1);
    ONSstat(5,y) = prctile(ONSy(:,y), 2.5, 1);
    ONSstat(6,y) = prctile(ONSy(:,y), 97.5, 1);
end

% Difference score analysis
TTPd = [TTP(:,3) - TTP(:,1), TTP(:,2) - TTP(:,1), TTP(:,3) - TTP(:,2), TTP(:,6) - TTP(:,4), TTP(:,5) - TTP(:,4), TTP(:,6) - TTP(:,5), TTP(:,9) - TTP(:,7), TTP(:,8) - TTP(:,7), TTP(:,9) - TTP(:,8)];
AMPd = [AMP(:,3) - AMP(:,1), AMP(:,2) - AMP(:,1), AMP(:,3) - AMP(:,2), AMP(:,6) - AMP(:,4), AMP(:,5) - AMP(:,4), AMP(:,6) - AMP(:,5), AMP(:,9) - AMP(:,7), AMP(:,8) - AMP(:,7), AMP(:,9) - AMP(:,8)];
ONSd = [ONS(:,3) - ONS(:,1), ONS(:,2) - ONS(:,1), ONS(:,3) - ONS(:,2), ONS(:,6) - ONS(:,4), ONS(:,5) - ONS(:,4), ONS(:,6) - ONS(:,5), ONS(:,9) - ONS(:,7), ONS(:,8) - ONS(:,7), ONS(:,9) - ONS(:,8)];

subjmTTPd = nanmean(TTPd,2);
grpmTTPd = nanmean(TTPd,'all');
subjmAMPd = nanmean(AMPd,2);
grpmAMPd = nanmean(AMPd,'all');
subjmONSd = nanmean(ONSd,2);
grpmONSd = nanmean(ONSd,'all'); 

for z = 1:size(TTP,1)
    TTPdy(z,:) = TTPd(z,:) - subjmTTPd(z,:) + grpmTTPd;
    AMPdy(z,:) = AMPd(z,:) - subjmAMPd(z,:) + grpmAMPd;
    ONSdy(z,:) = ONSd(z,:) - subjmONSd(z,:) + grpmONSd;
end

for y = 1:size(TTP,2)
    TTPdstat(1,y) = nanmean(TTPd(:,y),1); %means
    AMPdstat(1,y) = nanmean(AMPd(:,y),1);
    ONSdstat(1,y) = nanmean(ONSd(:,y),1);
    
    TTPdstat(2,y) = nanstd(TTPd(:,y),1); %sd
    AMPdstat(2,y) = nanstd(AMPd(:,y),1);
    ONSdstat(2,y) = nanstd(ONSd(:,y),1);
    
    TTPdstat(3,y) = prctile(TTPd(:,y), 2.5, 1); %ci
    TTPdstat(4,y) = prctile(TTPd(:,y), 97.5, 1);
    AMPdstat(3,y) = prctile(AMPd(:,y), 2.5, 1);
    AMPdstat(4,y) = prctile(AMPd(:,y), 97.5, 1);
    ONSdstat(3,y) = prctile(ONSd(:,y), 2.5, 1);
    ONSdstat(4,y) = prctile(ONSd(:,y), 97.5, 1);
end

TTPeffect(1,1) = TTPdstat(1,1) / sqrt( (nanvar(TTP(:,1),1) + nanvar(TTP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,2) = TTPdstat(1,2) / sqrt( (nanvar(TTP(:,1),1) + nanvar(TTP(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,3) = TTPdstat(1,3) / sqrt( (nanvar(TTP(:,2),1) + nanvar(TTP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,4) = TTPdstat(1,4) / sqrt( (nanvar(TTP(:,4),1) + nanvar(TTP(:,6),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,5) = TTPdstat(1,5) / sqrt( (nanvar(TTP(:,4),1) + nanvar(TTP(:,5),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,6) = TTPdstat(1,6) / sqrt( (nanvar(TTP(:,5),1) + nanvar(TTP(:,6),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,7) = TTPdstat(1,7) / sqrt( (nanvar(TTP(:,7),1) + nanvar(TTP(:,9),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,8) = TTPdstat(1,8) / sqrt( (nanvar(TTP(:,7),1) + nanvar(TTP(:,8),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,9) = TTPdstat(1,9) / sqrt( (nanvar(TTP(:,8),1) + nanvar(TTP(:,9),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));

AMPeffect(1,1) = AMPdstat(1,1) / sqrt( (nanvar(AMP(:,1),1) + nanvar(AMP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,2) = AMPdstat(1,2) / sqrt( (nanvar(AMP(:,1),1) + nanvar(AMP(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,3) = AMPdstat(1,3) / sqrt( (nanvar(AMP(:,2),1) + nanvar(AMP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,4) = AMPdstat(1,4) / sqrt( (nanvar(AMP(:,4),1) + nanvar(AMP(:,6),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,5) = AMPdstat(1,5) / sqrt( (nanvar(AMP(:,4),1) + nanvar(AMP(:,5),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,6) = AMPdstat(1,6) / sqrt( (nanvar(AMP(:,5),1) + nanvar(AMP(:,6),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,7) = AMPdstat(1,7) / sqrt( (nanvar(AMP(:,7),1) + nanvar(AMP(:,9),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,8) = AMPdstat(1,8) / sqrt( (nanvar(AMP(:,7),1) + nanvar(AMP(:,8),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,9) = AMPdstat(1,9) / sqrt( (nanvar(AMP(:,8),1) + nanvar(AMP(:,9),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));

ONSeffect(1,1) = ONSdstat(1,1) / sqrt( (nanvar(ONS(:,1),1) + nanvar(ONS(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,2) = ONSdstat(1,2) / sqrt( (nanvar(ONS(:,1),1) + nanvar(ONS(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,3) = ONSdstat(1,3) / sqrt( (nanvar(ONS(:,2),1) + nanvar(ONS(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,4) = ONSdstat(1,4) / sqrt( (nanvar(ONS(:,4),1) + nanvar(ONS(:,6),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,5) = ONSdstat(1,5) / sqrt( (nanvar(ONS(:,4),1) + nanvar(ONS(:,5),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,6) = ONSdstat(1,6) / sqrt( (nanvar(ONS(:,5),1) + nanvar(ONS(:,6),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,7) = ONSdstat(1,7) / sqrt( (nanvar(ONS(:,7),1) + nanvar(ONS(:,9),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,8) = ONSdstat(1,8) / sqrt( (nanvar(ONS(:,7),1) + nanvar(ONS(:,8),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,9) = ONSdstat(1,9) / sqrt( (nanvar(ONS(:,8),1) + nanvar(ONS(:,9),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));

TTPdstat(5,:) = TTPeffect;
AMPdstat(5,:) = AMPeffect;
ONSdstat(5,:) = ONSeffect;

cd(outputdir)
writematrix(TTPstat, sprintf('%s', 'BootStp', '_', 'TimeToPeak', '_', 'durationManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'BootStp', '_', 'Amplitude', '_', 'durationManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'BootStp', '_', 'PeakOnset', '_', 'durationManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'BootStp', '_', 'TimeToPeakD', '_', 'durationManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'BootStp', '_', 'AmplitudeD', '_', 'durationManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'BootStp', '_', 'PeakOnsetD', '_', 'durationManipulation', '.csv'))
%%
TTP = TTPstattableI; %table2array(stattable{x}(:,1:3));
AMP = AMPstattableI; %table2array(stattable{x}(:,4:6));
ONS = ONSstattableI; %table2array(stattable{x}(:,7:9));

subjmTTP = nanmean(TTP,2);
grpmTTP = nanmean(TTP,'all');
subjmAMP = nanmean(AMP,2);
grpmAMP = nanmean(AMP,'all');
subjmONS = nanmean(ONS,2);
grpmONS = nanmean(ONS,'all'); 

for z = 1:size(TTP,1)
    TTPy(z,:) = TTP(z,:) - subjmTTP(z,:) + grpmTTP;
    AMPy(z,:) = AMP(z,:) - subjmAMP(z,:) + grpmAMP;
    ONSy(z,:) = ONS(z,:) - subjmONS(z,:) + grpmONS;
end

for y = 1:size(TTP,2)
    TTPstat(1,y) = nanmean(TTP(:,y),1); %means
    AMPstat(1,y) = nanmean(AMP(:,y),1);
    ONSstat(1,y) = nanmean(ONS(:,y),1);
    
    TTPstat(2,y) = nanstd(TTP(:,y),1); %sd
    AMPstat(2,y) = nanstd(AMP(:,y),1);
    ONSstat(2,y) = nanstd(ONS(:,y),1);
    
    TTPstat(3,y) = prctile(TTP(:,y), 2.5, 1); %ci
    TTPstat(4,y) = prctile(TTP(:,y), 97.5, 1);
    AMPstat(3,y) = prctile(AMP(:,y), 2.5, 1);
    AMPstat(4,y) = prctile(AMP(:,y), 97.5, 1);
    ONSstat(3,y) = prctile(ONS(:,y), 2.5, 1);
    ONSstat(4,y) = prctile(ONS(:,y), 97.5, 1);

    TTPstat(5,y) = prctile(TTPy(:,y), 2.5, 1); %ci-within
    TTPstat(6,y) = prctile(TTPy(:,y), 97.5, 1);
    AMPstat(5,y) = prctile(AMPy(:,y), 2.5, 1);
    AMPstat(6,y) = prctile(AMPy(:,y), 97.5, 1);
    ONSstat(5,y) = prctile(ONSy(:,y), 2.5, 1);
    ONSstat(6,y) = prctile(ONSy(:,y), 97.5, 1);
end

% Difference score analysis
TTPd = [TTP(:,3) - TTP(:,1), TTP(:,2) - TTP(:,1), TTP(:,3) - TTP(:,2), TTP(:,6) - TTP(:,4), TTP(:,5) - TTP(:,4), TTP(:,6) - TTP(:,5), TTP(:,9) - TTP(:,7), TTP(:,8) - TTP(:,7), TTP(:,9) - TTP(:,8)];
AMPd = [AMP(:,3) - AMP(:,1), AMP(:,2) - AMP(:,1), AMP(:,3) - AMP(:,2), AMP(:,6) - AMP(:,4), AMP(:,5) - AMP(:,4), AMP(:,6) - AMP(:,5), AMP(:,9) - AMP(:,7), AMP(:,8) - AMP(:,7), AMP(:,9) - AMP(:,8)];
ONSd = [ONS(:,3) - ONS(:,1), ONS(:,2) - ONS(:,1), ONS(:,3) - ONS(:,2), ONS(:,6) - ONS(:,4), ONS(:,5) - ONS(:,4), ONS(:,6) - ONS(:,5), ONS(:,9) - ONS(:,7), ONS(:,8) - ONS(:,7), ONS(:,9) - ONS(:,8)];

subjmTTPd = nanmean(TTPd,2);
grpmTTPd = nanmean(TTPd,'all');
subjmAMPd = nanmean(AMPd,2);
grpmAMPd = nanmean(AMPd,'all');
subjmONSd = nanmean(ONSd,2);
grpmONSd = nanmean(ONSd,'all'); 

for z = 1:size(TTP,1)
    TTPdy(z,:) = TTPd(z,:) - subjmTTPd(z,:) + grpmTTPd;
    AMPdy(z,:) = AMPd(z,:) - subjmAMPd(z,:) + grpmAMPd;
    ONSdy(z,:) = ONSd(z,:) - subjmONSd(z,:) + grpmONSd;
end

for y = 1:size(TTP,2)
    TTPdstat(1,y) = nanmean(TTPd(:,y),1); %means
    AMPdstat(1,y) = nanmean(AMPd(:,y),1);
    ONSdstat(1,y) = nanmean(ONSd(:,y),1);
    
    TTPdstat(2,y) = nanstd(TTPd(:,y),1); %sd
    AMPdstat(2,y) = nanstd(AMPd(:,y),1);
    ONSdstat(2,y) = nanstd(ONSd(:,y),1);
    
    TTPdstat(3,y) = prctile(TTPd(:,y), 2.5, 1); %ci
    TTPdstat(4,y) = prctile(TTPd(:,y), 97.5, 1);
    AMPdstat(3,y) = prctile(AMPd(:,y), 2.5, 1);
    AMPdstat(4,y) = prctile(AMPd(:,y), 97.5, 1);
    ONSdstat(3,y) = prctile(ONSd(:,y), 2.5, 1);
    ONSdstat(4,y) = prctile(ONSd(:,y), 97.5, 1);
end

TTPeffect(1,1) = TTPdstat(1,1) / sqrt( (nanvar(TTP(:,1),1) + nanvar(TTP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,2) = TTPdstat(1,2) / sqrt( (nanvar(TTP(:,1),1) + nanvar(TTP(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,3) = TTPdstat(1,3) / sqrt( (nanvar(TTP(:,2),1) + nanvar(TTP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,4) = TTPdstat(1,4) / sqrt( (nanvar(TTP(:,4),1) + nanvar(TTP(:,6),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,5) = TTPdstat(1,5) / sqrt( (nanvar(TTP(:,4),1) + nanvar(TTP(:,5),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,6) = TTPdstat(1,6) / sqrt( (nanvar(TTP(:,5),1) + nanvar(TTP(:,6),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,7) = TTPdstat(1,7) / sqrt( (nanvar(TTP(:,7),1) + nanvar(TTP(:,9),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,8) = TTPdstat(1,8) / sqrt( (nanvar(TTP(:,7),1) + nanvar(TTP(:,8),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
TTPeffect(1,9) = TTPdstat(1,9) / sqrt( (nanvar(TTP(:,8),1) + nanvar(TTP(:,9),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));

AMPeffect(1,1) = AMPdstat(1,1) / sqrt( (nanvar(AMP(:,1),1) + nanvar(AMP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,2) = AMPdstat(1,2) / sqrt( (nanvar(AMP(:,1),1) + nanvar(AMP(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,3) = AMPdstat(1,3) / sqrt( (nanvar(AMP(:,2),1) + nanvar(AMP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,4) = AMPdstat(1,4) / sqrt( (nanvar(AMP(:,4),1) + nanvar(AMP(:,6),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,5) = AMPdstat(1,5) / sqrt( (nanvar(AMP(:,4),1) + nanvar(AMP(:,5),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,6) = AMPdstat(1,6) / sqrt( (nanvar(AMP(:,5),1) + nanvar(AMP(:,6),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,7) = AMPdstat(1,7) / sqrt( (nanvar(AMP(:,7),1) + nanvar(AMP(:,9),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,8) = AMPdstat(1,8) / sqrt( (nanvar(AMP(:,7),1) + nanvar(AMP(:,8),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));
AMPeffect(1,9) = AMPdstat(1,9) / sqrt( (nanvar(AMP(:,8),1) + nanvar(AMP(:,9),1)) / 2) * ((size(AMPd,1)-3) / (size(AMPd,1)-2.25)) * (sqrt((size(AMPd,1)-2) / size(AMPd,1)));

ONSeffect(1,1) = ONSdstat(1,1) / sqrt( (nanvar(ONS(:,1),1) + nanvar(ONS(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,2) = ONSdstat(1,2) / sqrt( (nanvar(ONS(:,1),1) + nanvar(ONS(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,3) = ONSdstat(1,3) / sqrt( (nanvar(ONS(:,2),1) + nanvar(ONS(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,4) = ONSdstat(1,4) / sqrt( (nanvar(ONS(:,4),1) + nanvar(ONS(:,6),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,5) = ONSdstat(1,5) / sqrt( (nanvar(ONS(:,4),1) + nanvar(ONS(:,5),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,6) = ONSdstat(1,6) / sqrt( (nanvar(ONS(:,5),1) + nanvar(ONS(:,6),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,7) = ONSdstat(1,7) / sqrt( (nanvar(ONS(:,7),1) + nanvar(ONS(:,9),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,8) = ONSdstat(1,8) / sqrt( (nanvar(ONS(:,7),1) + nanvar(ONS(:,8),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));
ONSeffect(1,9) = ONSdstat(1,9) / sqrt( (nanvar(ONS(:,8),1) + nanvar(ONS(:,9),1)) / 2) * ((size(ONSd,1)-3) / (size(ONSd,1)-2.25)) * (sqrt((size(ONSd,1)-2) / size(ONSd,1)));

TTPdstat(5,:) = TTPeffect;
AMPdstat(5,:) = AMPeffect;
ONSdstat(5,:) = ONSeffect;

cd(outputdir)
writematrix(TTPstat, sprintf('%s', 'BootStp', '_', 'TimeToPeak', '_', 'intensityManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'BootStp', '_', 'Amplitude', '_', 'intensityManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'BootStp', '_', 'PeakOnset', '_', 'intensityManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'BootStp', '_', 'TimeToPeakD', '_', 'intensityManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'BootStp', '_', 'AmplitudeD', '_', 'intensityManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'BootStp', '_', 'PeakOnsetD', '_', 'intensityManipulation', '.csv'))
