clearvars()
cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\out\vis_newscript'); % Path to data matrices output by PerSubjProc
addpath 'C:\Users\Alvin\iCloudDrive\tron_analysis'
eventtype = "V"; % Set ROI type (A or V)
filedir = dir('V_F*.mat'); % Set filename search criteria
varNames = {'TTP1', 'TTP2', 'TTP3', 'AMP1', 'AMP2', 'AMP3', 'ONS1', 'ONS2', 'ONS3', 'TTPERA1', 'TTPERA2', 'TTPERA3', 'AMPERA1', 'AMPERA2', 'AMPERA3',};
varTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
stattable{1} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{2} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{3} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);

stattable{4} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{5} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{6} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);

for fileidx = 1:numel(filedir)
    clearvars('-except','varNames', 'varTypes', 'stattable',... 
        'eventtype', 'filedir', 'customvolsperevent', 'volsperevent', 'fileidx', 'AllStats', 'SubjFittedHRFs', 'SubjFittedFullHRFs', 'SubjERAs');
    filename = sprintf('%s', [filedir(fileidx).folder, '/', filedir(fileidx).name]);
    load(filename);
    subjname = [filedir(fileidx).name(4:10)]; % Data import and subject name definition

    conddef = {["Length 0.1 Int 0.1" "Length 0.1 Int 0.3" "Length 0.1 Int 0.9";... 
                "Length 0.3 Int 0.1" "Length 0.3 Int 0.3" "Length 0.3 Int 0.9";... 
                "Length 0.9 Int 0.1" "Length 0.9 Int 0.3" "Length 0.9 Int 0.9"];...
                ["Length 0.1 Int 0.01" "Length 0.1 Int 0.1" "Length 0.1 Int 1";... 
                "Length 0.3 Int 0.01" "Length 0.3 Int 0.1" "Length 0.3 Int 1";...
                "Length 0.9 Int 0.01" "Length 0.9 Int 0.1" "Length 0.9 Int 1"]};
            
    for lengths = 1:3
        for ints = 1:3
            DelayMatrix{lengths, ints} = FittingMatrix{lengths, ints}.delay;
            AmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.amplitude;
            DispersionMatrix{lengths, ints} = FittingMatrix{lengths, ints}.dispersion;
            OnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.onset;
            PeakTimeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxIndexInSeconds; % Best peak delay estimate
            PeakAmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxValue; % Best peak amplitude estimate
            FWHMMatrix{lengths, ints} = FittingMatrix{lengths, ints}.FWHMInSeconds; % Best dispersion estimate
            PeakOnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.OnsetTimeInSeconds; % Best peak onset estimate
        end
    end
            
    stattable{1}(fileidx,1:3) = PeakTimeMatrix(:,1)';
    stattable{1}(fileidx,4:6) = PeakAmplitudeMatrix(:,1)';
    stattable{1}(fileidx, 7:9) = PeakOnsetMatrix(:,1)';
    
    stattable{2}(fileidx,1:3) = PeakTimeMatrix(:,2)';
    stattable{2}(fileidx,4:6) = PeakAmplitudeMatrix(:,2)';
    stattable{2}(fileidx, 7:9) = PeakOnsetMatrix(:,2)';

    stattable{3}(fileidx,1:3) = PeakTimeMatrix(:,3)';
    stattable{3}(fileidx,4:6) = PeakAmplitudeMatrix(:,3)';
    stattable{3}(fileidx, 7:9) = PeakOnsetMatrix(:,3)';

    stattable{4}(fileidx,1:3) = PeakTimeMatrix(1,:);
    stattable{4}(fileidx,4:6) = PeakAmplitudeMatrix(1,:);
    stattable{4}(fileidx, 7:9) = PeakOnsetMatrix(1,:);

    
    stattable{5}(fileidx,1:3) = PeakTimeMatrix(2,:);
    stattable{5}(fileidx,4:6) = PeakAmplitudeMatrix(2,:);
    stattable{5}(fileidx, 7:9) = PeakOnsetMatrix(2,:);

    stattable{6}(fileidx,1:3) = PeakTimeMatrix(3,:);
    stattable{6}(fileidx,4:6) = PeakAmplitudeMatrix(3,:);
    stattable{6}(fileidx, 7:9) = PeakOnsetMatrix(3,:);

end

%%
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
    
    TTPstat(2,y) = nanstd(TTPy(:,y),1); %sd
    AMPstat(2,y) = nanstd(AMPy(:,y),1);
    ONSstat(2,y) = nanstd(ONSy(:,y),1);

    TTPstat(3,y) = sqrt(nanvar(TTPy(:,y),1) * 1/size(TTPy,1)); %sem
    AMPstat(3,y) = sqrt(nanvar(AMPy(:,y),1) * 1/size(AMPy,1));
    ONSstat(3,y) = sqrt(nanvar(ONSy(:,y),1) * 1/size(ONSy,1));

    TTPstat(4,y) = nanmean(TTP(:,y),1) - (1.96 * sqrt(nanvar(TTPy(:,y),1) * 1/size(TTP,1))); %ci
    TTPstat(5,y) = nanmean(TTP(:,y),1) + (1.96 * sqrt(nanvar(TTPy(:,y),1) * 1/size(TTP,1)));
    AMPstat(4,y) = nanmean(AMP(:,y),1) - (1.96 * sqrt(nanvar(AMPy(:,y),1) * 1/size(AMP,1))); 
    AMPstat(5,y) = nanmean(AMP(:,y),1) + (1.96 * sqrt(nanvar(AMPy(:,y),1) * 1/size(AMP,1)));
    ONSstat(4,y) = nanmean(ONS(:,y),1) - (1.96 * sqrt(nanvar(ONSy(:,y),1) * 1/size(ONS,1)));
    ONSstat(5,y) = nanmean(ONS(:,y),1) + (1.96 * sqrt(nanvar(ONSy(:,y),1) * 1/size(ONS,1)));
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
    TTPdstat(2,y) = nanstd(TTPdy(:,y),1); %sd
    AMPdstat(2,y) = nanstd(AMPdy(:,y),1);
    ONSdstat(2,y) = nanstd(ONSdy(:,y),1);
    TTPdstat(3,y) = sqrt(nanvar(TTPdy(:,y),1) * 1/size(TTPd,1)); %sem
    AMPdstat(3,y) = sqrt(nanvar(AMPdy(:,y),1) * 1/size(AMPd,1));
    ONSdstat(3,y) = sqrt(nanvar(ONSdy(:,y),1) * 1/size(ONSd,1));
    TTPdstat(4,y) = nanmean(TTPd(:,y),1) - (1.96 * sqrt(nanvar(TTPdy(:,y),1) * 1/size(TTPd,1))); %ci
    TTPdstat(5,y) = nanmean(TTPd(:,y),1) + (1.96 * sqrt(nanvar(TTPdy(:,y),1) * 1/size(TTPd,1)));
    AMPdstat(4,y) = nanmean(AMPd(:,y),1) - (1.96 * sqrt(nanvar(AMPdy(:,y),1) * 1/size(AMPd,1))); 
    AMPdstat(5,y) = nanmean(AMPd(:,y),1) + (1.96 * sqrt(nanvar(AMPdy(:,y),1) * 1/size(AMPd,1)));
    ONSdstat(4,y) = nanmean(ONSd(:,y),1) - (1.96 * sqrt(nanvar(ONSdy(:,y),1) * 1/size(ONSd,1)));
    ONSdstat(5,y) = nanmean(ONSd(:,y),1) + (1.96 * sqrt(nanvar(ONSdy(:,y),1) * 1/size(ONSd,1)));
    
    
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

TTPdstat(6,:) = TTPeffect;
AMPdstat(6,:) = AMPeffect;
ONSdstat(6,:) = ONSeffect;

writematrix(TTPstat, sprintf('%s', 'TimeToPeak', '_', 'durationManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'Amplitude', '_', 'durationManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'PeakOnset', '_', 'durationManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'TimeToPeakD', '_', 'durationManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'AmplitudeD', '_', 'durationManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'PeakOnsetD', '_', 'durationManipulation', '.csv'))


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

    TTPstat(2,y) = nanstd(TTPy(:,y),1); %sd
    AMPstat(2,y) = nanstd(AMPy(:,y),1);
    ONSstat(2,y) = nanstd(ONSy(:,y),1);

    TTPstat(3,y) = sqrt(nanvar(TTPy(:,y),1) * 1/size(TTPy,1)); %sem
    AMPstat(3,y) = sqrt(nanvar(AMPy(:,y),1) * 1/size(AMPy,1));
    ONSstat(3,y) = sqrt(nanvar(ONSy(:,y),1) * 1/size(ONSy,1));

    TTPstat(4,y) = nanmean(TTP(:,y),1) - (1.96 * sqrt(nanvar(TTPy(:,y),1) * 1/size(TTP,1))); %ci
    TTPstat(5,y) = nanmean(TTP(:,y),1) + (1.96 * sqrt(nanvar(TTPy(:,y),1) * 1/size(TTP,1)));
    AMPstat(4,y) = nanmean(AMP(:,y),1) - (1.96 * sqrt(nanvar(AMPy(:,y),1) * 1/size(AMP,1))); 
    AMPstat(5,y) = nanmean(AMP(:,y),1) + (1.96 * sqrt(nanvar(AMPy(:,y),1) * 1/size(AMP,1)));
    ONSstat(4,y) = nanmean(ONS(:,y),1) - (1.96 * sqrt(nanvar(ONSy(:,y),1) * 1/size(ONS,1)));
    ONSstat(5,y) = nanmean(ONS(:,y),1) + (1.96 * sqrt(nanvar(ONSy(:,y),1) * 1/size(ONS,1)));
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
    TTPdstat(2,y) = nanstd(TTPdy(:,y),1); %sd
    AMPdstat(2,y) = nanstd(AMPdy(:,y),1);
    ONSdstat(2,y) = nanstd(ONSdy(:,y),1);
    TTPdstat(3,y) = sqrt(nanvar(TTPdy(:,y),1) * 1/size(TTPd,1)); %sem
    AMPdstat(3,y) = sqrt(nanvar(AMPdy(:,y),1) * 1/size(AMPd,1));
    ONSdstat(3,y) = sqrt(nanvar(ONSdy(:,y),1) * 1/size(ONSd,1));
    TTPdstat(4,y) = nanmean(TTPd(:,y),1) - (1.96 * sqrt(nanvar(TTPdy(:,y),1) * 1/size(TTPd,1))); %ci
    TTPdstat(5,y) = nanmean(TTPd(:,y),1) + (1.96 * sqrt(nanvar(TTPdy(:,y),1) * 1/size(TTPd,1)));
    AMPdstat(4,y) = nanmean(AMPd(:,y),1) - (1.96 * sqrt(nanvar(AMPdy(:,y),1) * 1/size(AMPd,1))); 
    AMPdstat(5,y) = nanmean(AMPd(:,y),1) + (1.96 * sqrt(nanvar(AMPdy(:,y),1) * 1/size(AMPd,1)));
    ONSdstat(4,y) = nanmean(ONSd(:,y),1) - (1.96 * sqrt(nanvar(ONSdy(:,y),1) * 1/size(ONSd,1)));
    ONSdstat(5,y) = nanmean(ONSd(:,y),1) + (1.96 * sqrt(nanvar(ONSdy(:,y),1) * 1/size(ONSd,1)));
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

TTPdstat(6,:) = TTPeffect;
AMPdstat(6,:) = AMPeffect;
ONSdstat(6,:) = ONSeffect;

writematrix(TTPstat, sprintf('%s', 'TimeToPeak', '_', 'intensityManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'Amplitude', '_', 'intensityManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'PeakOnset', '_', 'intensityManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'TimeToPeakD', '_', 'intensityManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'AmplitudeD', '_', 'intensityManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'PeakOnsetD', '_', 'intensityManipulation', '.csv'))

%% DURATION EFFECTS WITHOUT WITHIN-SUBJ CORRECTIONS
TTP = TTPstattableL; %table2array(stattable{x}(:,1:3));
AMP = AMPstattableL; %table2array(stattable{x}(:,4:6));
ONS = ONSstattableL; %table2array(stattable{x}(:,7:9));

for y = 1:size(TTP,2)
    TTPstat(1,y) = nanmean(TTP(:,y),1); %means
    AMPstat(1,y) = nanmean(AMP(:,y),1);
    ONSstat(1,y) = nanmean(ONS(:,y),1);
    
    TTPstat(2,y) = nanstd(TTP(:,y),1); %sd
    AMPstat(2,y) = nanstd(AMP(:,y),1);
    ONSstat(2,y) = nanstd(ONS(:,y),1);

    TTPstat(3,y) = sqrt(nanvar(TTP(:,y),1) * 1/size(TTP,1)); %sem
    AMPstat(3,y) = sqrt(nanvar(AMP(:,y),1) * 1/size(AMP,1));
    ONSstat(3,y) = sqrt(nanvar(ONS(:,y),1) * 1/size(ONS,1));

    TTPstat(4,y) = nanmean(TTP(:,y),1) - (1.96 * sqrt(nanvar(TTP(:,y),1) * 1/size(TTP,1))); %ci
    TTPstat(5,y) = nanmean(TTP(:,y),1) + (1.96 * sqrt(nanvar(TTP(:,y),1) * 1/size(TTP,1)));
    AMPstat(4,y) = nanmean(AMP(:,y),1) - (1.96 * sqrt(nanvar(AMP(:,y),1) * 1/size(AMP,1))); 
    AMPstat(5,y) = nanmean(AMP(:,y),1) + (1.96 * sqrt(nanvar(AMP(:,y),1) * 1/size(AMP,1)));
    ONSstat(4,y) = nanmean(ONS(:,y),1) - (1.96 * sqrt(nanvar(ONS(:,y),1) * 1/size(ONS,1)));
    ONSstat(5,y) = nanmean(ONS(:,y),1) + (1.96 * sqrt(nanvar(ONS(:,y),1) * 1/size(ONS,1)));

end

% Difference score analysis
TTPd = [TTP(:,3) - TTP(:,1), TTP(:,2) - TTP(:,1), TTP(:,3) - TTP(:,2), TTP(:,6) - TTP(:,4), TTP(:,5) - TTP(:,4), TTP(:,6) - TTP(:,5), TTP(:,9) - TTP(:,7), TTP(:,8) - TTP(:,7), TTP(:,9) - TTP(:,8)];
AMPd = [AMP(:,3) - AMP(:,1), AMP(:,2) - AMP(:,1), AMP(:,3) - AMP(:,2), AMP(:,6) - AMP(:,4), AMP(:,5) - AMP(:,4), AMP(:,6) - AMP(:,5), AMP(:,9) - AMP(:,7), AMP(:,8) - AMP(:,7), AMP(:,9) - AMP(:,8)];
ONSd = [ONS(:,3) - ONS(:,1), ONS(:,2) - ONS(:,1), ONS(:,3) - ONS(:,2), ONS(:,6) - ONS(:,4), ONS(:,5) - ONS(:,4), ONS(:,6) - ONS(:,5), ONS(:,9) - ONS(:,7), ONS(:,8) - ONS(:,7), ONS(:,9) - ONS(:,8)];

for y = 1:size(TTP,2)
    TTPdstat(1,y) = nanmean(TTPd(:,y),1); %means
    AMPdstat(1,y) = nanmean(AMPd(:,y),1);
    ONSdstat(1,y) = nanmean(ONSd(:,y),1);
    TTPdstat(2,y) = nanstd(TTPd(:,y),1); %sd
    AMPdstat(2,y) = nanstd(AMPd(:,y),1);
    ONSdstat(2,y) = nanstd(ONSd(:,y),1);
    TTPdstat(3,y) = sqrt(nanvar(TTPd(:,y),1) * 1/size(TTPd,1)); %sem
    AMPdstat(3,y) = sqrt(nanvar(AMPd(:,y),1) * 1/size(AMPd,1));
    ONSdstat(3,y) = sqrt(nanvar(ONSd(:,y),1) * 1/size(ONSd,1));
    TTPdstat(4,y) = nanmean(TTPd(:,y),1) - (1.96 * sqrt(nanvar(TTPd(:,y),1) * 1/size(TTPd,1))); %ci
    TTPdstat(5,y) = nanmean(TTPd(:,y),1) + (1.96 * sqrt(nanvar(TTPd(:,y),1) * 1/size(TTPd,1)));
    AMPdstat(4,y) = nanmean(AMPd(:,y),1) - (1.96 * sqrt(nanvar(AMPd(:,y),1) * 1/size(AMPd,1))); 
    AMPdstat(5,y) = nanmean(AMPd(:,y),1) + (1.96 * sqrt(nanvar(AMPd(:,y),1) * 1/size(AMPd,1)));
    ONSdstat(4,y) = nanmean(ONSd(:,y),1) - (1.96 * sqrt(nanvar(ONSd(:,y),1) * 1/size(ONSd,1)));
    ONSdstat(5,y) = nanmean(ONSd(:,y),1) + (1.96 * sqrt(nanvar(ONSd(:,y),1) * 1/size(ONSd,1)));
    
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

TTPdstat(6,:) = TTPeffect;
AMPdstat(6,:) = AMPeffect;
ONSdstat(6,:) = ONSeffect;

writematrix(TTPstat, sprintf('%s', 'NoWithinCorrections', '_', 'TimeToPeak', '_', 'durationManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'NoWithinCorrections', '_', 'Amplitude', '_', 'durationManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'NoWithinCorrections', '_', 'PeakOnset', '_', 'durationManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'NoWithinCorrections', '_', 'TimeToPeakD', '_', 'durationManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'NoWithinCorrections', '_', 'AmplitudeD', '_', 'durationManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'NoWithinCorrections', '_', 'PeakOnsetD', '_', 'durationManipulation', '.csv'))

%% INTENSITY EFFECTS WITHOUT WITHIN-SUBJ CORRECTIONS
TTP = TTPstattableI; %table2array(stattable{x}(:,1:3));
AMP = AMPstattableI; %table2array(stattable{x}(:,4:6));
ONS = ONSstattableI; %table2array(stattable{x}(:,7:9));


for y = 1:size(TTP,2)
    TTPstat(1,y) = nanmean(TTP(:,y),1); %means
    AMPstat(1,y) = nanmean(AMP(:,y),1);
    ONSstat(1,y) = nanmean(ONS(:,y),1);
    
    TTPstat(2,y) = nanstd(TTP(:,y),1); %sd
    AMPstat(2,y) = nanstd(AMP(:,y),1);
    ONSstat(2,y) = nanstd(ONS(:,y),1);
    
    TTPstat(3,y) = sqrt(nanvar(TTP(:,y),1) * 1/size(TTP,1)); %sem
    AMPstat(3,y) = sqrt(nanvar(AMP(:,y),1) * 1/size(AMP,1));
    ONSstat(3,y) = sqrt(nanvar(ONS(:,y),1) * 1/size(ONS,1));

    TTPstat(4,y) = nanmean(TTP(:,y),1) - (1.96 * sqrt(nanvar(TTP(:,y),1) * 1/size(TTP,1))); %ci
    TTPstat(5,y) = nanmean(TTP(:,y),1) + (1.96 * sqrt(nanvar(TTP(:,y),1) * 1/size(TTP,1)));
    AMPstat(4,y) = nanmean(AMP(:,y),1) - (1.96 * sqrt(nanvar(AMP(:,y),1) * 1/size(AMP,1))); 
    AMPstat(5,y) = nanmean(AMP(:,y),1) + (1.96 * sqrt(nanvar(AMP(:,y),1) * 1/size(AMP,1)));
    ONSstat(4,y) = nanmean(ONS(:,y),1) - (1.96 * sqrt(nanvar(ONS(:,y),1) * 1/size(ONS,1)));
    ONSstat(5,y) = nanmean(ONS(:,y),1) + (1.96 * sqrt(nanvar(ONS(:,y),1) * 1/size(ONS,1)));

end

% Difference score analysis
TTPd = [TTP(:,3) - TTP(:,1), TTP(:,2) - TTP(:,1), TTP(:,3) - TTP(:,2), TTP(:,6) - TTP(:,4), TTP(:,5) - TTP(:,4), TTP(:,6) - TTP(:,5), TTP(:,9) - TTP(:,7), TTP(:,8) - TTP(:,7), TTP(:,9) - TTP(:,8)];
AMPd = [AMP(:,3) - AMP(:,1), AMP(:,2) - AMP(:,1), AMP(:,3) - AMP(:,2), AMP(:,6) - AMP(:,4), AMP(:,5) - AMP(:,4), AMP(:,6) - AMP(:,5), AMP(:,9) - AMP(:,7), AMP(:,8) - AMP(:,7), AMP(:,9) - AMP(:,8)];
ONSd = [ONS(:,3) - ONS(:,1), ONS(:,2) - ONS(:,1), ONS(:,3) - ONS(:,2), ONS(:,6) - ONS(:,4), ONS(:,5) - ONS(:,4), ONS(:,6) - ONS(:,5), ONS(:,9) - ONS(:,7), ONS(:,8) - ONS(:,7), ONS(:,9) - ONS(:,8)];

for y = 1:size(TTP,2)
    TTPdstat(1,y) = nanmean(TTPd(:,y),1); %means
    AMPdstat(1,y) = nanmean(AMPd(:,y),1);
    ONSdstat(1,y) = nanmean(ONSd(:,y),1);
    TTPdstat(2,y) = nanstd(TTPd(:,y),1); %sd
    AMPdstat(2,y) = nanstd(AMPd(:,y),1);
    ONSdstat(2,y) = nanstd(ONSd(:,y),1);
    TTPdstat(3,y) = sqrt(nanvar(TTPd(:,y),1) * 1/size(TTPd,1)); %sem
    AMPdstat(3,y) = sqrt(nanvar(AMPd(:,y),1) * 1/size(AMPd,1));
    ONSdstat(3,y) = sqrt(nanvar(ONSd(:,y),1) * 1/size(ONSd,1));
    TTPdstat(4,y) = nanmean(TTPd(:,y),1) - (1.96 * sqrt(nanvar(TTPd(:,y),1) * 1/size(TTPd,1))); %ci
    TTPdstat(5,y) = nanmean(TTPd(:,y),1) + (1.96 * sqrt(nanvar(TTPd(:,y),1) * 1/size(TTPd,1)));
    AMPdstat(4,y) = nanmean(AMPd(:,y),1) - (1.96 * sqrt(nanvar(AMPd(:,y),1) * 1/size(AMPd,1))); 
    AMPdstat(5,y) = nanmean(AMPd(:,y),1) + (1.96 * sqrt(nanvar(AMPd(:,y),1) * 1/size(AMPd,1)));
    ONSdstat(4,y) = nanmean(ONSd(:,y),1) - (1.96 * sqrt(nanvar(ONSd(:,y),1) * 1/size(ONSd,1)));
    ONSdstat(5,y) = nanmean(ONSd(:,y),1) + (1.96 * sqrt(nanvar(ONSd(:,y),1) * 1/size(ONSd,1)));
    
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

TTPdstat(6,:) = TTPeffect;
AMPdstat(6,:) = AMPeffect;
ONSdstat(6,:) = ONSeffect;

writematrix(TTPstat, sprintf('%s', 'NoWithinCorrections', '_', 'TimeToPeak', '_', 'intensityManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'NoWithinCorrections', '_', 'Amplitude', '_', 'intensityManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s','NoWithinCorrections', '_',  'PeakOnset', '_', 'intensityManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s','NoWithinCorrections', '_',  'TimeToPeakD', '_', 'intensityManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'NoWithinCorrections', '_', 'AmplitudeD', '_', 'intensityManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'NoWithinCorrections', '_', 'PeakOnsetD', '_', 'intensityManipulation', '.csv'))
