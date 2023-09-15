clearvars();

% Per-subject ROI timecourse event-related averaging, linear gamma function
% fitting and point estimate computation

% Written Sep 2018 by Alvin Wong for TRoN Project

cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\Data_Pub'); % Path to data matrices output by PerSubjProc
addpath 'C:\Users\Alvin\iCloudDrive\tron_analysis'
eventtype = "V"; % Set ROI type (A or V)
filedir = dir('V1*.mat'); % Set filename search criteria

varNames = {'TTP1', 'TTP2', 'TTP3', 'AMP1', 'AMP2', 'AMP3', 'ONS1', 'ONS2', 'ONS3', 'TTPERA1', 'TTPERA2', 'TTPERA3', 'AMPERA1', 'AMPERA2', 'AMPERA3',};
varTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
stattable{1} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{2} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{3} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);

stattable{4} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{5} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{6} = table('Size', [numel(filedir), 3*5], 'VariableTypes', varTypes, 'VariableNames', varNames);

 conddef = {["Length 0.1 Int 0.1" "Length 0.1 Int 0.3" "Length 0.1 Int 0.9";... 
                "Length 0.3 Int 0.1" "Length 0.3 Int 0.3" "Length 0.3 Int 0.9";... 
                "Length 0.9 Int 0.1" "Length 0.9 Int 0.3" "Length 0.9 Int 0.9"];...
                ["Length 0.1 Int 0.01" "Length 0.1 Int 0.1" "Length 0.1 Int 1";... 
                "Length 0.3 Int 0.01" "Length 0.3 Int 0.1" "Length 0.3 Int 1";...
                "Length 0.9 Int 0.01" "Length 0.9 Int 0.1" "Length 0.9 Int 1"]};


        
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
    if eventtype == 'V'
        conddefinition = conddef{2};
        eventlengths = [0.1 0.3 0.9];
        eventintensities = [0.01 0.1 1];
        roiname = 'V1';
        outputdir = "C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\out\vis_newscript";
    elseif eventtype == 'A'
        conddefinition = conddef{1};
        eventlengths = [0.1 0.3 0.9];
        eventintensities = [0.1 0.3 1];
        roiname = 'A1';
        outputdir = "C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\out\aud_newscript";
    end            

    for conC = 1:numel(ProtoMatrix) % baseline computation
        BaselineVector{conC} = ProtoMatrix{conC}(13, :); % event onsets at index 10
        for i = 1:size(ProtoMatrix{conC}, 2)
            ProtoMatrix{conC}(:, i) = ((ProtoMatrix{conC}(:, i) - BaselineVector{conC}(1,i)) / BaselineVector{conC}(1,i)) * 100;
        end
        EventAvg{conC} = nanmean(ProtoMatrix{conC}(10:34,:), 2);
        PeakERA{conC} = max(EventAvg{conC});
        TimeToPeakERA{conC} = find(abs(EventAvg{conC}-PeakERA{conC}) < 0.001) * 0.625;
        TimeToPeakERA{conC} = TimeToPeakERA{conC}(1);
        ErrMargin{conC} = (nanstd(ProtoMatrix{conC}(10:34,:),0,2)/sqrt(size(ProtoMatrix{conC},2)))*1.96;
    end

        EventAvg = reshape(EventAvg, [3 3]);
        ErrMargin = reshape(ErrMargin, [3 3]);
        PeakERA = reshape(PeakERA, [3 3]);
        TimeToPeakERA = reshape(TimeToPeakERA, [3 3]);

    % Per-condition linear gamma function fitting for each subject
    for lengths = 1:size(ProtoMatrix,1)
        for ints = 1:size(ProtoMatrix,2)
            [FittingMatrix{lengths, ints}] = GammaFit(ProtoMatrix{lengths, ints}(10:34,:)', 0.625, 15);
            DelayMatrix{lengths, ints} = FittingMatrix{lengths, ints}.delay;
            AmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.amplitude;
            DispersionMatrix{lengths, ints} = FittingMatrix{lengths, ints}.dispersion;
            OnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.onset;
            PeakTimeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxIndexInSeconds; % Best peak delay estimate
            PeakAmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxValue; % Best peak amplitude estimate
            FWHMMatrix{lengths, ints} = FittingMatrix{lengths, ints}.FWHMInSeconds; % Best dispersion estimate
            PeakOnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.OnsetTimeInSeconds; % Best peak onset estimate
            
           FittedHRFs{lengths,ints} = FittingMatrix{lengths,ints}.hrf;
           FittedFullHRFs{lengths,ints} = FittingMatrix{lengths,ints}.fullhrf;
        end
    end
    
    cd(outputdir)
    save(sprintf('%s', [eventtype, '_', subjname, '_FittingMatrix', '.mat']),'FittingMatrix');

    AllStats{fileidx,1} = [{PeakTimeMatrix}, {PeakAmplitudeMatrix}, {PeakERA}, {TimeToPeakERA}, {FWHMMatrix}, {PeakOnsetMatrix}];
    SubjFittedHRFs{fileidx} = FittedHRFs;
    SubjFittedFullHRFs{fileidx} = FittedFullHRFs;
    SubjERAs{fileidx} = EventAvg;
    
    stattable{1}(fileidx,1:3) = PeakTimeMatrix(:,1)';
    stattable{1}(fileidx,4:6) = PeakAmplitudeMatrix(:,1)';
    stattable{1}(fileidx, 7:9) = PeakOnsetMatrix(:,1)';
    stattable{1}(fileidx, 10:12) = TimeToPeakERA(:,1)';
    stattable{1}(fileidx, 13:15) = PeakERA(:,1)';
    
    stattable{2}(fileidx,1:3) = PeakTimeMatrix(:,2)';
    stattable{2}(fileidx,4:6) = PeakAmplitudeMatrix(:,2)';
    stattable{2}(fileidx, 7:9) = PeakOnsetMatrix(:,2)';
    stattable{2}(fileidx, 10:12) = TimeToPeakERA(:,2)';
    stattable{2}(fileidx, 13:15) = PeakERA(:,2)';
    
    stattable{3}(fileidx,1:3) = PeakTimeMatrix(:,3)';
    stattable{3}(fileidx,4:6) = PeakAmplitudeMatrix(:,3)';
    stattable{3}(fileidx, 7:9) = PeakOnsetMatrix(:,3)';
    stattable{3}(fileidx, 10:12) = TimeToPeakERA(:,3)';
    stattable{3}(fileidx, 13:15) = PeakERA(:,3)';
    
    stattable{4}(fileidx,1:3) = PeakTimeMatrix(1,:);
    stattable{4}(fileidx,4:6) = PeakAmplitudeMatrix(1,:);
    stattable{4}(fileidx, 7:9) = PeakOnsetMatrix(1,:);
    stattable{4}(fileidx, 10:12) = TimeToPeakERA(1,:);
    stattable{4}(fileidx, 13:15) = PeakERA(1,:);
    
    stattable{5}(fileidx,1:3) = PeakTimeMatrix(2,:);
    stattable{5}(fileidx,4:6) = PeakAmplitudeMatrix(2,:);
    stattable{5}(fileidx, 7:9) = PeakOnsetMatrix(2,:);
    stattable{5}(fileidx, 10:12) = TimeToPeakERA(2,:);
    stattable{5}(fileidx, 13:15) = PeakERA(2,:);
    
    stattable{6}(fileidx,1:3) = PeakTimeMatrix(3,:);
    stattable{6}(fileidx,4:6) = PeakAmplitudeMatrix(3,:);
    stattable{6}(fileidx, 7:9) = PeakOnsetMatrix(3,:);
    stattable{6}(fileidx, 10:12) = TimeToPeakERA(3,:);
    stattable{6}(fileidx, 13:15) = PeakERA(3,:);

end
%%
for subj = 1:numel(SubjFittedHRFs)
    for lengths=1:3
        for ints=1:3
            HRFmatrix{lengths,ints}(subj,:) = SubjFittedHRFs{subj}{lengths,ints};
            FullHRFmatrix{lengths,ints}(subj,:) = SubjFittedFullHRFs{subj}{lengths,ints};
            ERAmatrix{lengths,ints}{subj,:} = SubjERAs{subj}{lengths,ints}';
        end
    end
end

for lengths=1:3
    for ints=1:3
        HRFmatrix{lengths,ints} = nanmean(HRFmatrix{lengths, ints},1);
        FullHRFmatrix{lengths,ints} = nanmean(FullHRFmatrix{lengths, ints},1);
        ERAmatrix{lengths,ints} = nanmean(cell2mat(ERAmatrix{lengths,ints}),1);
    end
end

TTPstattableL = [table2array(stattable{1}(:,1:3)), table2array(stattable{2}(:,1:3)), table2array(stattable{3}(:,1:3))];
TTPstattableI = [table2array(stattable{4}(:,1:3)), table2array(stattable{5}(:,1:3)), table2array(stattable{6}(:,1:3))];
AMPstattableL = [table2array(stattable{1}(:,4:6)), table2array(stattable{2}(:,4:6)), table2array(stattable{3}(:,4:6))];
AMPstattableI = [table2array(stattable{4}(:,4:6)), table2array(stattable{5}(:,4:6)), table2array(stattable{6}(:,4:6))];
ONSstattableL = [table2array(stattable{1}(:,7:9)), table2array(stattable{2}(:,7:9)), table2array(stattable{3}(:,7:9))];
ONSstattableI = [table2array(stattable{4}(:,7:9)), table2array(stattable{5}(:,7:9)), table2array(stattable{6}(:,7:9))];
TTPERAstattableL = [table2array(stattable{1}(:,10:12)), table2array(stattable{2}(:,10:12)), table2array(stattable{3}(:,10:12))];
TTPERAstattableI = [table2array(stattable{4}(:,10:12)), table2array(stattable{5}(:,10:12)), table2array(stattable{6}(:,10:12))];
ERAstattableL = [table2array(stattable{1}(:,13:15)), table2array(stattable{2}(:,13:15)), table2array(stattable{3}(:,13:15))];
ERAstattableI = [table2array(stattable{4}(:,13:15)), table2array(stattable{5}(:,13:15)), table2array(stattable{6}(:,13:15))];

comparisons = ["durations", "durations", "durations", "intensities", "intensities", "intensities"];
levels = ["intensity1", "intensity2", "intensity3", "duration1", "duration2", "duration3"];


%% DURATION EFFECTS WITH WITHIN-SUBJ CORRECTIONS

TTP = TTPstattableL; %table2array(stattable{x}(:,1:3));
AMP = AMPstattableL; %table2array(stattable{x}(:,4:6));
ONS = ONSstattableL; %table2array(stattable{x}(:,7:9));
TTPERA = TTPERAstattableL;
AMPERA = ERAstattableL;

subjmTTP = nanmean(TTP,2);
grpmTTP = nanmean(TTP,'all');
subjmAMP = nanmean(AMP,2);
grpmAMP = nanmean(AMP,'all');
subjmONS = nanmean(ONS,2);
grpmONS = nanmean(ONS,'all'); 
subjmTTPERA = nanmean(TTPERA,2);
grpmTTPERA = nanmean(TTPERA,'all');
subjmAMPERA = nanmean(AMPERA,2);
grpmAMPERA = nanmean(AMPERA,'all');

for z = 1:size(TTP,1)
    TTPy(z,:) = TTP(z,:) - subjmTTP(z,:) + grpmTTP;
    AMPy(z,:) = AMP(z,:) - subjmAMP(z,:) + grpmAMP;
    ONSy(z,:) = ONS(z,:) - subjmONS(z,:) + grpmONS;
    TTPERAy(z,:) = TTPERA(z,:) - subjmTTPERA(z,:) + grpmTTPERA;
    AMPERAy(z,:) = AMPERA(z,:) - subjmAMPERA(z,:) + grpmAMPERA;
end

for y = 1:size(TTP,2)
    TTPstat(1,y) = nanmean(TTP(:,y),1); %means
    AMPstat(1,y) = nanmean(AMP(:,y),1);
    ONSstat(1,y) = nanmean(ONS(:,y),1);
    TTPERAstat(1,y) = nanmean(TTPERA(:,y),1); 
    AMPERAstat(1,y) = nanmean(AMPERA(:,y),1);
    
    TTPstat(2,y) = nanstd(TTPy(:,y),1); %sd
    AMPstat(2,y) = nanstd(AMPy(:,y),1);
    ONSstat(2,y) = nanstd(ONSy(:,y),1);
    TTPERAstat(2,y) = nanstd(TTPERAy(:,y),1);
    AMPERAstat(2,y) = nanstd(AMPERAy(:,y),1);
    
    TTPstat(3,y) = sqrt(nanvar(TTPy(:,y),1) * 1/size(TTPy,1)); %sem
    AMPstat(3,y) = sqrt(nanvar(AMPy(:,y),1) * 1/size(AMPy,1));
    ONSstat(3,y) = sqrt(nanvar(ONSy(:,y),1) * 1/size(ONSy,1));
    TTPERAstat(3,y) = sqrt(nanvar(TTPERAy(:,y),1) * 1/size(TTPERAy,1));
    AMPERAstat(3,y) = sqrt(nanvar(AMPERAy(:,y),1) * 1/size(AMPERAy,1));
    
    TTPstat(4,y) = nanmean(TTP(:,y),1) - (1.96 * sqrt(nanvar(TTPy(:,y),1) * 1/size(TTP,1))); %ci
    TTPstat(5,y) = nanmean(TTP(:,y),1) + (1.96 * sqrt(nanvar(TTPy(:,y),1) * 1/size(TTP,1)));
    AMPstat(4,y) = nanmean(AMP(:,y),1) - (1.96 * sqrt(nanvar(AMPy(:,y),1) * 1/size(AMP,1))); 
    AMPstat(5,y) = nanmean(AMP(:,y),1) + (1.96 * sqrt(nanvar(AMPy(:,y),1) * 1/size(AMP,1)));
    ONSstat(4,y) = nanmean(ONS(:,y),1) - (1.96 * sqrt(nanvar(ONSy(:,y),1) * 1/size(ONS,1)));
    ONSstat(5,y) = nanmean(ONS(:,y),1) + (1.96 * sqrt(nanvar(ONSy(:,y),1) * 1/size(ONS,1)));
    TTPERAstat(4,y) = nanmean(TTPERA(:,y),1) - (1.96 * sqrt(nanvar(TTPERAy(:,y),1) * 1/size(TTPERA,1)));
    TTPERAstat(5,y) = nanmean(TTPERA(:,y),1) + (1.96 * sqrt(nanvar(TTPERAy(:,y),1) * 1/size(TTPERA,1)));
    AMPERAstat(4,y) = nanmean(AMPERA(:,y),1) - (1.96 * sqrt(nanvar(AMPERAy(:,y),1) * 1/size(AMPERA,1))); 
    AMPERAstat(5,y) = nanmean(AMPERA(:,y),1) + (1.96 * sqrt(nanvar(AMPERAy(:,y),1) * 1/size(AMPERA,1)));
end

% Difference score analysis
TTPd = [TTP(:,3) - TTP(:,1), TTP(:,2) - TTP(:,1), TTP(:,3) - TTP(:,2), TTP(:,6) - TTP(:,4), TTP(:,5) - TTP(:,4), TTP(:,6) - TTP(:,5), TTP(:,9) - TTP(:,7), TTP(:,8) - TTP(:,7), TTP(:,9) - TTP(:,8)];
AMPd = [AMP(:,3) - AMP(:,1), AMP(:,2) - AMP(:,1), AMP(:,3) - AMP(:,2), AMP(:,6) - AMP(:,4), AMP(:,5) - AMP(:,4), AMP(:,6) - AMP(:,5), AMP(:,9) - AMP(:,7), AMP(:,8) - AMP(:,7), AMP(:,9) - AMP(:,8)];
ONSd = [ONS(:,3) - ONS(:,1), ONS(:,2) - ONS(:,1), ONS(:,3) - ONS(:,2), ONS(:,6) - ONS(:,4), ONS(:,5) - ONS(:,4), ONS(:,6) - ONS(:,5), ONS(:,9) - ONS(:,7), ONS(:,8) - ONS(:,7), ONS(:,9) - ONS(:,8)];
TTPERAd = [TTPERA(:,3) - TTPERA(:,1), TTPERA(:,2) - TTPERA(:,1), TTPERA(:,3) - TTPERA(:,2), TTPERA(:,6) - TTPERA(:,4), TTPERA(:,5) - TTPERA(:,4), TTPERA(:,6) - TTPERA(:,5), TTPERA(:,9) - TTPERA(:,7), TTPERA(:,8) - TTPERA(:,7), TTPERA(:,9) - TTPERA(:,8)];
AMPERAd = [AMPERA(:,3) - AMPERA(:,1), AMPERA(:,2) - AMPERA(:,1), AMPERA(:,3) - AMPERA(:,2), AMPERA(:,6) - AMPERA(:,4), AMPERA(:,5) - AMPERA(:,4), AMPERA(:,6) - AMPERA(:,5), AMPERA(:,9) - AMPERA(:,7), AMPERA(:,8) - AMPERA(:,7), AMPERA(:,9) - AMPERA(:,8)];


subjmTTPd = nanmean(TTPd,2);
grpmTTPd = nanmean(TTPd,'all');
subjmAMPd = nanmean(AMPd,2);
grpmAMPd = nanmean(AMPd,'all');
subjmONSd = nanmean(ONSd,2);
grpmONSd = nanmean(ONSd,'all'); 

subjmTTPERAd = nanmean(TTPERAd,2);
grpmTTPERAd = nanmean(TTPERAd,'all');
subjmAMPERAd = nanmean(AMPERAd,2);
grpmAMPERAd = nanmean(AMPERAd,'all');

for z = 1:size(TTP,1)
    TTPdy(z,:) = TTPd(z,:) - subjmTTPd(z,:) + grpmTTPd;
    AMPdy(z,:) = AMPd(z,:) - subjmAMPd(z,:) + grpmAMPd;
    ONSdy(z,:) = ONSd(z,:) - subjmONSd(z,:) + grpmONSd;
    TTPERAdy(z,:) = TTPERAd(z,:) - subjmTTPERAd(z,:) + grpmTTPERAd;
    AMPERAdy(z,:) = AMPERAd(z,:) - subjmAMPERAd(z,:) + grpmAMPERAd;
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
    
    TTPERAdstat(1,y) = nanmean(TTPERAd(:,y),1); %means
    AMPERAdstat(1,y) = nanmean(AMPERAd(:,y),1);
    
    TTPERAdstat(2,y) = nanstd(TTPERAdy(:,y),1); %sd
    AMPERAdstat(2,y) = nanstd(AMPERAdy(:,y),1);
    
    TTPERAdstat(3,y) = sqrt(nanvar(TTPERAdy(:,y),1) * 1/size(TTPERAd,1)); %sem
    AMPERAdstat(3,y) = sqrt(nanvar(AMPERAdy(:,y),1) * 1/size(AMPERAd,1));
   
    TTPERAdstat(4,y) = nanmean(TTPERAd(:,y),1) - (1.96 * sqrt(nanvar(TTPERAdy(:,y),1) * 1/size(TTPERAd,1))); %ci
    TTPERAdstat(5,y) = nanmean(TTPERAd(:,y),1) + (1.96 * sqrt(nanvar(TTPERAdy(:,y),1) * 1/size(TTPERAd,1)));
    AMPERAdstat(4,y) = nanmean(AMPERAd(:,y),1) - (1.96 * sqrt(nanvar(AMPERAdy(:,y),1) * 1/size(AMPERAd,1))); 
    AMPERAdstat(5,y) = nanmean(AMPERAd(:,y),1) + (1.96 * sqrt(nanvar(AMPERAdy(:,y),1) * 1/size(AMPERAd,1)));
    
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

TTPERAeffect(1,1) = TTPERAdstat(1,1) / sqrt( (nanvar(TTPERA(:,1),1) + nanvar(TTPERA(:,3),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,2) = TTPERAdstat(1,2) / sqrt( (nanvar(TTPERA(:,1),1) + nanvar(TTPERA(:,2),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,3) = TTPERAdstat(1,3) / sqrt( (nanvar(TTPERA(:,2),1) + nanvar(TTPERA(:,3),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,4) = TTPERAdstat(1,4) / sqrt( (nanvar(TTPERA(:,4),1) + nanvar(TTPERA(:,6),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,5) = TTPERAdstat(1,5) / sqrt( (nanvar(TTPERA(:,4),1) + nanvar(TTPERA(:,5),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,6) = TTPERAdstat(1,6) / sqrt( (nanvar(TTPERA(:,5),1) + nanvar(TTPERA(:,6),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,7) = TTPERAdstat(1,7) / sqrt( (nanvar(TTPERA(:,7),1) + nanvar(TTPERA(:,9),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,8) = TTPERAdstat(1,8) / sqrt( (nanvar(TTPERA(:,7),1) + nanvar(TTPERA(:,8),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,9) = TTPERAdstat(1,9) / sqrt( (nanvar(TTPERA(:,8),1) + nanvar(TTPERA(:,9),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));

AMPERAeffect(1,1) = AMPERAdstat(1,1) / sqrt( (nanvar(AMPERA(:,1),1) + nanvar(AMPERA(:,3),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
AMPERAeffect(1,2) = AMPERAdstat(1,2) / sqrt( (nanvar(AMPERA(:,1),1) + nanvar(AMPERA(:,2),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
AMPERAeffect(1,3) = AMPERAdstat(1,3) / sqrt( (nanvar(AMPERA(:,2),1) + nanvar(AMPERA(:,3),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
AMPERAeffect(1,4) = AMPERAdstat(1,4) / sqrt( (nanvar(AMPERA(:,4),1) + nanvar(AMPERA(:,6),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,5) = AMPERAdstat(1,5) / sqrt( (nanvar(AMPERA(:,4),1) + nanvar(AMPERA(:,5),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,6) = AMPERAdstat(1,6) / sqrt( (nanvar(AMPERA(:,5),1) + nanvar(AMPERA(:,6),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,7) = AMPERAdstat(1,7) / sqrt( (nanvar(AMPERA(:,7),1) + nanvar(AMPERA(:,9),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,8) = AMPERAdstat(1,8) / sqrt( (nanvar(AMPERA(:,7),1) + nanvar(AMPERA(:,8),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,9) = AMPERAdstat(1,9) / sqrt( (nanvar(AMPERA(:,8),1) + nanvar(AMPERA(:,9),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));


TTPdstat(6,:) = TTPeffect;
AMPdstat(6,:) = AMPeffect;
ONSdstat(6,:) = ONSeffect;
TTPERAdstat(6,:) = TTPERAeffect;
AMPERAdstat(6,:) = AMPERAeffect;

writematrix(TTPstat, sprintf('%s', 'TimeToPeak', '_', 'durationManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'Amplitude', '_', 'durationManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'PeakOnset', '_', 'durationManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'TimeToPeakD', '_', 'durationManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'AmplitudeD', '_', 'durationManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'PeakOnsetD', '_', 'durationManipulation', '.csv'))

writematrix(TTPERAstat, sprintf('%s', 'TimeToPeakERA', '_', 'durationManipulation', '.csv'))
writematrix(AMPERAstat, sprintf('%s', 'PeakERA', '_', 'durationManipulation', '.csv'))
writematrix(TTPERAdstat, sprintf('%s', 'TimeToPeakERAD', '_', 'durationManipulation', '.csv'))
writematrix(AMPERAdstat, sprintf('%s', 'PeakERAD', '_', 'durationManipulation', '.csv'))

%% INTENSITY EFFECTS WITH WITHIN-SUBJ CORRECTIONS
TTP = TTPstattableI; %table2array(stattable{x}(:,1:3));
AMP = AMPstattableI; %table2array(stattable{x}(:,4:6));
ONS = ONSstattableI; %table2array(stattable{x}(:,7:9));
TTPERA = TTPERAstattableI;
AMPERA = ERAstattableI;

subjmTTP = nanmean(TTP,2);
grpmTTP = nanmean(TTP,'all');
subjmAMP = nanmean(AMP,2);
grpmAMP = nanmean(AMP,'all');
subjmONS = nanmean(ONS,2);
grpmONS = nanmean(ONS,'all'); 
subjmTTPERA = nanmean(TTPERA,2);
grpmTTPERA = nanmean(TTPERA,'all');
subjmAMPERA = nanmean(AMPERA,2);
grpmAMPERA = nanmean(AMPERA,'all');

for z = 1:size(TTP,1)
    TTPy(z,:) = TTP(z,:) - subjmTTP(z,:) + grpmTTP;
    AMPy(z,:) = AMP(z,:) - subjmAMP(z,:) + grpmAMP;
    ONSy(z,:) = ONS(z,:) - subjmONS(z,:) + grpmONS;
    TTPERAy(z,:) = TTPERA(z,:) - subjmTTPERA(z,:) + grpmTTPERA;
    AMPERAy(z,:) = AMPERA(z,:) - subjmAMPERA(z,:) + grpmAMPERA;
end

for y = 1:size(TTP,2)
    TTPstat(1,y) = nanmean(TTP(:,y),1); %means
    AMPstat(1,y) = nanmean(AMP(:,y),1);
    ONSstat(1,y) = nanmean(ONS(:,y),1);
    TTPERAstat(1,y) = nanmean(TTPERA(:,y),1); 
    AMPERAstat(1,y) = nanmean(AMPERA(:,y),1);
    
    TTPstat(2,y) = nanstd(TTPy(:,y),1); %sd
    AMPstat(2,y) = nanstd(AMPy(:,y),1);
    ONSstat(2,y) = nanstd(ONSy(:,y),1);
    TTPERAstat(2,y) = nanstd(TTPERAy(:,y),1);
    AMPERAstat(2,y) = nanstd(AMPERAy(:,y),1);
    
    TTPstat(3,y) = sqrt(nanvar(TTPy(:,y),1) * 1/size(TTPy,1)); %sem
    AMPstat(3,y) = sqrt(nanvar(AMPy(:,y),1) * 1/size(AMPy,1));
    ONSstat(3,y) = sqrt(nanvar(ONSy(:,y),1) * 1/size(ONSy,1));
    TTPERAstat(3,y) = sqrt(nanvar(TTPERAy(:,y),1) * 1/size(TTPERAy,1));
    AMPERAstat(3,y) = sqrt(nanvar(AMPERAy(:,y),1) * 1/size(AMPERAy,1));
    
    TTPstat(4,y) = nanmean(TTP(:,y),1) - (1.96 * sqrt(nanvar(TTPy(:,y),1) * 1/size(TTP,1))); %ci
    TTPstat(5,y) = nanmean(TTP(:,y),1) + (1.96 * sqrt(nanvar(TTPy(:,y),1) * 1/size(TTP,1)));
    AMPstat(4,y) = nanmean(AMP(:,y),1) - (1.96 * sqrt(nanvar(AMPy(:,y),1) * 1/size(AMP,1))); 
    AMPstat(5,y) = nanmean(AMP(:,y),1) + (1.96 * sqrt(nanvar(AMPy(:,y),1) * 1/size(AMP,1)));
    ONSstat(4,y) = nanmean(ONS(:,y),1) - (1.96 * sqrt(nanvar(ONSy(:,y),1) * 1/size(ONS,1)));
    ONSstat(5,y) = nanmean(ONS(:,y),1) + (1.96 * sqrt(nanvar(ONSy(:,y),1) * 1/size(ONS,1)));
    TTPERAstat(4,y) = nanmean(TTPERA(:,y),1) - (1.96 * sqrt(nanvar(TTPERAy(:,y),1) * 1/size(TTPERA,1)));
    TTPERAstat(5,y) = nanmean(TTPERA(:,y),1) + (1.96 * sqrt(nanvar(TTPERAy(:,y),1) * 1/size(TTPERA,1)));
    AMPERAstat(4,y) = nanmean(AMPERA(:,y),1) - (1.96 * sqrt(nanvar(AMPERAy(:,y),1) * 1/size(AMPERA,1))); 
    AMPERAstat(5,y) = nanmean(AMPERA(:,y),1) + (1.96 * sqrt(nanvar(AMPERAy(:,y),1) * 1/size(AMPERA,1)));
end

% Difference score analysis
TTPd = [TTP(:,3) - TTP(:,1), TTP(:,2) - TTP(:,1), TTP(:,3) - TTP(:,2), TTP(:,6) - TTP(:,4), TTP(:,5) - TTP(:,4), TTP(:,6) - TTP(:,5), TTP(:,9) - TTP(:,7), TTP(:,8) - TTP(:,7), TTP(:,9) - TTP(:,8)];
AMPd = [AMP(:,3) - AMP(:,1), AMP(:,2) - AMP(:,1), AMP(:,3) - AMP(:,2), AMP(:,6) - AMP(:,4), AMP(:,5) - AMP(:,4), AMP(:,6) - AMP(:,5), AMP(:,9) - AMP(:,7), AMP(:,8) - AMP(:,7), AMP(:,9) - AMP(:,8)];
ONSd = [ONS(:,3) - ONS(:,1), ONS(:,2) - ONS(:,1), ONS(:,3) - ONS(:,2), ONS(:,6) - ONS(:,4), ONS(:,5) - ONS(:,4), ONS(:,6) - ONS(:,5), ONS(:,9) - ONS(:,7), ONS(:,8) - ONS(:,7), ONS(:,9) - ONS(:,8)];
TTPERAd = [TTPERA(:,3) - TTPERA(:,1), TTPERA(:,2) - TTPERA(:,1), TTPERA(:,3) - TTPERA(:,2), TTPERA(:,6) - TTPERA(:,4), TTPERA(:,5) - TTPERA(:,4), TTPERA(:,6) - TTPERA(:,5), TTPERA(:,9) - TTPERA(:,7), TTPERA(:,8) - TTPERA(:,7), TTPERA(:,9) - TTPERA(:,8)];
AMPERAd = [AMPERA(:,3) - AMPERA(:,1), AMPERA(:,2) - AMPERA(:,1), AMPERA(:,3) - AMPERA(:,2), AMPERA(:,6) - AMPERA(:,4), AMPERA(:,5) - AMPERA(:,4), AMPERA(:,6) - AMPERA(:,5), AMPERA(:,9) - AMPERA(:,7), AMPERA(:,8) - AMPERA(:,7), AMPERA(:,9) - AMPERA(:,8)];


subjmTTPd = nanmean(TTPd,2);
grpmTTPd = nanmean(TTPd,'all');
subjmAMPd = nanmean(AMPd,2);
grpmAMPd = nanmean(AMPd,'all');
subjmONSd = nanmean(ONSd,2);
grpmONSd = nanmean(ONSd,'all'); 

subjmTTPERAd = nanmean(TTPERAd,2);
grpmTTPERAd = nanmean(TTPERAd,'all');
subjmAMPERAd = nanmean(AMPERAd,2);
grpmAMPERAd = nanmean(AMPERAd,'all');

for z = 1:size(TTP,1)
    TTPdy(z,:) = TTPd(z,:) - subjmTTPd(z,:) + grpmTTPd;
    AMPdy(z,:) = AMPd(z,:) - subjmAMPd(z,:) + grpmAMPd;
    ONSdy(z,:) = ONSd(z,:) - subjmONSd(z,:) + grpmONSd;
    TTPERAdy(z,:) = TTPERAd(z,:) - subjmTTPERAd(z,:) + grpmTTPERAd;
    AMPERAdy(z,:) = AMPERAd(z,:) - subjmAMPERAd(z,:) + grpmAMPERAd;
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
    
    TTPERAdstat(1,y) = nanmean(TTPERAd(:,y),1); %means
    AMPERAdstat(1,y) = nanmean(AMPERAd(:,y),1);
    
    TTPERAdstat(2,y) = nanstd(TTPERAdy(:,y),1); %sd
    AMPERAdstat(2,y) = nanstd(AMPERAdy(:,y),1);
    
    TTPERAdstat(3,y) = sqrt(nanvar(TTPERAdy(:,y),1) * 1/size(TTPERAd,1)); %sem
    AMPERAdstat(3,y) = sqrt(nanvar(AMPERAdy(:,y),1) * 1/size(AMPERAd,1));
   
    TTPERAdstat(4,y) = nanmean(TTPERAd(:,y),1) - (1.96 * sqrt(nanvar(TTPERAdy(:,y),1) * 1/size(TTPERAd,1))); %ci
    TTPERAdstat(5,y) = nanmean(TTPERAd(:,y),1) + (1.96 * sqrt(nanvar(TTPERAdy(:,y),1) * 1/size(TTPERAd,1)));
    AMPERAdstat(4,y) = nanmean(AMPERAd(:,y),1) - (1.96 * sqrt(nanvar(AMPERAdy(:,y),1) * 1/size(AMPERAd,1))); 
    AMPERAdstat(5,y) = nanmean(AMPERAd(:,y),1) + (1.96 * sqrt(nanvar(AMPERAdy(:,y),1) * 1/size(AMPERAd,1)));
    
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

TTPERAeffect(1,1) = TTPERAdstat(1,1) / sqrt( (nanvar(TTPERA(:,1),1) + nanvar(TTPERA(:,3),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,2) = TTPERAdstat(1,2) / sqrt( (nanvar(TTPERA(:,1),1) + nanvar(TTPERA(:,2),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,3) = TTPERAdstat(1,3) / sqrt( (nanvar(TTPERA(:,2),1) + nanvar(TTPERA(:,3),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,4) = TTPERAdstat(1,4) / sqrt( (nanvar(TTPERA(:,4),1) + nanvar(TTPERA(:,6),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,5) = TTPERAdstat(1,5) / sqrt( (nanvar(TTPERA(:,4),1) + nanvar(TTPERA(:,5),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,6) = TTPERAdstat(1,6) / sqrt( (nanvar(TTPERA(:,5),1) + nanvar(TTPERA(:,6),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,7) = TTPERAdstat(1,7) / sqrt( (nanvar(TTPERA(:,7),1) + nanvar(TTPERA(:,9),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,8) = TTPERAdstat(1,8) / sqrt( (nanvar(TTPERA(:,7),1) + nanvar(TTPERA(:,8),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
TTPERAeffect(1,9) = TTPERAdstat(1,9) / sqrt( (nanvar(TTPERA(:,8),1) + nanvar(TTPERA(:,9),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));

AMPERAeffect(1,1) = AMPERAdstat(1,1) / sqrt( (nanvar(AMPERA(:,1),1) + nanvar(AMPERA(:,3),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
AMPERAeffect(1,2) = AMPERAdstat(1,2) / sqrt( (nanvar(AMPERA(:,1),1) + nanvar(AMPERA(:,2),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
AMPERAeffect(1,3) = AMPERAdstat(1,3) / sqrt( (nanvar(AMPERA(:,2),1) + nanvar(AMPERA(:,3),1)) / 2) * ((size(TTPERAd,1)-3) / (size(TTPERAd,1)-2.25)) * (sqrt((size(TTPERAd,1)-2) / size(TTPERAd,1)));
AMPERAeffect(1,4) = AMPERAdstat(1,4) / sqrt( (nanvar(AMPERA(:,4),1) + nanvar(AMPERA(:,6),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,5) = AMPERAdstat(1,5) / sqrt( (nanvar(AMPERA(:,4),1) + nanvar(AMPERA(:,5),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,6) = AMPERAdstat(1,6) / sqrt( (nanvar(AMPERA(:,5),1) + nanvar(AMPERA(:,6),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,7) = AMPERAdstat(1,7) / sqrt( (nanvar(AMPERA(:,7),1) + nanvar(AMPERA(:,9),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,8) = AMPERAdstat(1,8) / sqrt( (nanvar(AMPERA(:,7),1) + nanvar(AMPERA(:,8),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));
AMPERAeffect(1,9) = AMPERAdstat(1,9) / sqrt( (nanvar(AMPERA(:,8),1) + nanvar(AMPERA(:,9),1)) / 2) * ((size(AMPERAd,1)-3) / (size(AMPERAd,1)-2.25)) * (sqrt((size(AMPERAd,1)-2) / size(AMPERAd,1)));

TTPdstat(6,:) = TTPeffect;
AMPdstat(6,:) = AMPeffect;
ONSdstat(6,:) = ONSeffect;
TTPERAdstat(6,:) = TTPERAeffect;
AMPERAdstat(6,:) = AMPERAeffect;

writematrix(TTPstat, sprintf('%s', 'TimeToPeak', '_', 'intensityManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'Amplitude', '_', 'intensityManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'PeakOnset', '_', 'intensityManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'TimeToPeakD', '_', 'intensityManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'AmplitudeD', '_', 'intensityManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'PeakOnsetD', '_', 'intensityManipulation', '.csv'))

writematrix(TTPERAstat, sprintf('%s', 'TimeToPeakERA', '_', 'intensityManipulation', '.csv'))
writematrix(AMPERAstat, sprintf('%s', 'PeakERA', '_', 'intensityManipulation', '.csv'))
writematrix(TTPERAdstat, sprintf('%s', 'TimeToPeakERAD', '_', 'intensityManipulation', '.csv'))
writematrix(AMPERAdstat, sprintf('%s', 'PeakERAD', '_', 'intensityManipulation', '.csv'))


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
