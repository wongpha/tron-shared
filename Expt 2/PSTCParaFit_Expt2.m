clearvars();

% Per-subject ROI timecourse event-related averaging, linear gamma function
% fitting and point estimate computation

% Written Sep 2018 by Alvin Wong for TRoN Project

cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\PSubjDataT'); % Path to data matrices output by PerSubjProc
addpath 'C:\Users\alvinw\iCloudDrive\Documents\tron_analysis'
eventtype = "V"; % Set ROI type (A or V)
filedir = dir('V1*.mat'); % Set filename search criteria
varNames = {'TTP1', 'TTP2', 'TTP3', 'AMP1', 'AMP2', 'AMP3', 'ONS1', 'ONS2', 'ONS3','TTPa', 'AMPa', 'ONSa'};
varTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
stattable = table('Size', [numel(filedir), 3*4], 'VariableTypes', varTypes, 'VariableNames', varNames);
for fileidx = 1:numel(filedir)
    clearvars('-except', 'varNames', 'varTypes', 'stattable', 'eventtype', 'filedir', 'customvolsperevent', 'volsperevent', 'fileidx', 'AllStats', 'SubjFittedHRFs', 'SubjFittedFullHRFs', 'SubjERAs');
    filename = sprintf('%s', [filedir(fileidx).folder, '/', filedir(fileidx).name]);
    load(filename);
    subjname = [filedir(fileidx).name(4:10)]; % Data import and subject name definition

    conddef = {["NA" "Length 1.2 Int 0.3" "NA"; "Length 1.0 Int 0.9" "Length 1.2 Int 0.9" "Length 1.8 Int 0.9"]';...
    ["NA" "Length 1.2 Int 0.1" "NA"; "Length 1.0 Int 1" "Length 1.2 Int 1" "Length 1.8 Int 1"]'};

                
    if eventtype == 'V' % V1 timecourse processing to build protomatrix
        conddefinition = conddef{2};
        eventlengths = [1.0 1.2 1.8];
        eventintensities = [0.1 1];
        roiname = 'V1';
        
        for datacnt = 1:numel(CurrentV)
        VOnset(:, datacnt) = CurrentV(datacnt).VisualEvs(1:end);
        VLength(:,datacnt) = CurrentV(datacnt).VisualLengths;
        VInt(:,datacnt) = CurrentV(datacnt).VisualInts;
        end
        
        for datacnt = 1:numel(CurrentV)
        for evcnt = 1:size(VOnset(:,datacnt),1)
            VTC(:,evcnt,datacnt) = CurrentV(datacnt).RawData((VOnset(evcnt,datacnt)-4):(VOnset(evcnt,datacnt)+15));
        end
        end

        for datacnt = 1:size(VLength,2)
        for evcnt = 1:size(VLength(:,datacnt),1)
            for Lcon = 1:numel(eventlengths)
                Lsearch{Lcon}(evcnt, datacnt) = eq(VLength(evcnt, datacnt),eventlengths(Lcon));
            end
        end
        for evcnt = 1:size(VInt(:,datacnt),1)
            for Icon = 1:numel(eventintensities)
                Isearch{Icon}(evcnt, datacnt) = eq(VInt(evcnt, datacnt),eventintensities(Icon));
            end
        end
        end

        ProtoMatrix{3, 2} = [];

        for Lcon = 1:numel(eventlengths)
        for Icon = 1:numel(eventintensities)
            ConditionMap{Lcon, Icon} = [eventlengths(Lcon) eventintensities(Icon)];
            for datacnt = 1:size(Lsearch{Lcon},2)
                matcnt = 1;
                for evcnt = 1:size(Lsearch{Lcon}(:,datacnt),1)
                    if Lsearch{Lcon}(evcnt, datacnt) == 1 && Isearch{Icon}(evcnt, datacnt) == 1
                        ProtoMatrix{Lcon, Icon}(1:size(VTC,1), matcnt, datacnt) = VTC(1:end, evcnt,datacnt);
                        matcnt = matcnt+1;
                    end
                end
            end
        end
        end

    elseif eventtype == 'A' % A1 timecourse processing to build protomatrix
        conddefinition = conddef{1};
        eventlengths = [1.0 1.2 1.8];
        eventintensities = [0.3 1];

        roiname = 'A1';
        
        for datacnt = 1:numel(CurrentA)
        AOnset(:, datacnt) = CurrentA(datacnt).AuditoryEvs(1:end);
        ALength(:,datacnt) = CurrentA(datacnt).AuditoryLengths;
        AInt(:,datacnt) = CurrentA(datacnt).AuditoryInts;
        end
        
        for datacnt = 1:numel(CurrentA)
        for evcnt = 1:size(AOnset(:,datacnt),1)
            ATC(:,evcnt,datacnt) = CurrentA(datacnt).RawData((AOnset(evcnt,datacnt)-4):(AOnset(evcnt,datacnt)+15));
        end
        end

        for datacnt = 1:size(ALength,2)
        for evcnt = 1:size(ALength(:,datacnt),1)
            for Lcon = 1:numel(eventlengths)
                Lsearch{Lcon}(evcnt, datacnt) = eq(ALength(evcnt, datacnt),eventlengths(Lcon));
            end
        end
        for evcnt = 1:size(AInt(:,datacnt),1)
            for Icon = 1:numel(eventintensities)
                Isearch{Icon}(evcnt, datacnt) = eq(AInt(evcnt, datacnt),eventintensities(Icon));
            end
        end
        end

        ProtoMatrix{3, 2} = [];

        for Lcon = 1:numel(eventlengths)
        for Icon = 1:numel(eventintensities)
            ConditionMap{Lcon, Icon} = [eventlengths(Lcon) eventintensities(Icon)];
            for datacnt = 1:size(Lsearch{Lcon},2)
                matcnt = 1;
                for evcnt = 1:size(Lsearch{Lcon}(:,datacnt),1)
                    if Lsearch{Lcon}(evcnt, datacnt) == 1 && Isearch{Icon}(evcnt, datacnt) == 1
                        ProtoMatrix{Lcon, Icon}(1:size(ATC,1), matcnt, datacnt) = ATC(1:end, evcnt,datacnt);
                        matcnt = matcnt+1;
                    end
                end
            end
        end
        end

    end
    
    for conC = 1:numel(ProtoMatrix) % convert 3D cells in protomatrix to 2D cells
        ProtoMatrix{conC} = reshape(ProtoMatrix{conC}, [20, size(ProtoMatrix{conC},2)*size(ProtoMatrix{conC},3)]);
    end
    
    for conC = 1:numel(ProtoMatrix) % baseline computation
        BaselineVector{conC} = ProtoMatrix{conC}(4, :);
        for i = 1:size(ProtoMatrix{conC}, 2)
            ProtoMatrix{conC}(:, i) = ((ProtoMatrix{conC}(:, i) - BaselineVector{conC}(1,i)) / BaselineVector{conC}(1,i)) * 100;
        end
        EventAvg{conC} = nanmean(ProtoMatrix{conC}(5:20,:), 2);
        PeakERA{conC} = max(EventAvg{conC});
        TimeToPeakERA{conC} = find(abs(EventAvg{conC}-PeakERA{conC}) < 0.001) * 1;
        %TimeToPeakERA{conC} = TimeToPeakERA{conC}(1);
        ErrMargin{conC} = (nanstd(ProtoMatrix{conC}(5:20,:),0,2)/sqrt(size(ProtoMatrix{conC},2)))*1.96;
    end

        EventAvg = reshape(EventAvg, [3 2]);
        ErrMargin = reshape(ErrMargin, [3 2]);
        PeakERA = reshape(PeakERA, [3 2]);
        TimeToPeakERA = reshape(TimeToPeakERA, [3 2]);
    
    cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\PSubjData\Out'); % save all event averages to MATs
    save(sprintf('%s', [subjname, '_', roiname, '_ERAvgs']), 'ProtoMatrix', 'EventAvg', 'ErrMargin', 'PeakERA', 'TimeToPeakERA');
    
    % Per-condition linear gamma function fitting for each subject
    for lengths = 1:size(ProtoMatrix,1)
        for ints = 2 %1:size(ProtoMatrix,2)
            [FittingMatrix{lengths, ints}] = GammaFit(ProtoMatrix{lengths, ints}(5:20,:)', 1, 15);
            DelayMatrix{lengths, ints} = FittingMatrix{lengths, ints}.delay;
            AmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.amplitude;
            DispersionMatrix{lengths, ints} = FittingMatrix{lengths, ints}.dispersion;
            OnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.onset;
            PeakTimeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxIndexInSeconds; % Best peak delay estimate
            PeakAmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxValue; % Best peak amplitude estimate
            FWHMMatrix{lengths, ints} = FittingMatrix{lengths, ints}.FWHMInSeconds; % Best dispersion estimate
            PeakOnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.OnsetTimeInSeconds; % Best peak onset estimate
            
            hold on;
            
           FittedHRFs{lengths,ints} = FittingMatrix{lengths,ints}.hrf;
           FittedFullHRFs{lengths,ints} = FittingMatrix{lengths,ints}.fullhrf;
        end
    end
    for lengths = 2%:size(ProtoMatrix,1)
        for ints = 1 %1:size(ProtoMatrix,2)
            [FittingMatrix{lengths, ints}] = GammaFit(ProtoMatrix{lengths, ints}(5:20,:)', 1, 15);
            DelayMatrix{lengths, ints} = FittingMatrix{lengths, ints}.delay;
            AmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.amplitude;
            DispersionMatrix{lengths, ints} = FittingMatrix{lengths, ints}.dispersion;
            OnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.onset;
            PeakTimeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxIndexInSeconds; % Best peak delay estimate
            PeakAmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxValue; % Best peak amplitude estimate
            FWHMMatrix{lengths, ints} = FittingMatrix{lengths, ints}.FWHMInSeconds; % Best dispersion estimate
            PeakOnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.OnsetTimeInSeconds; % Best peak onset estimate
            
            hold on;
            
           FittedHRFs{lengths,ints} = FittingMatrix{lengths,ints}.hrf;
           FittedFullHRFs{lengths,ints} = FittingMatrix{lengths,ints}.fullhrf;
        end
    end
    
    save(sprintf('%s', [eventtype, '_', subjname, '_FittingMatrix', '.mat']),'FittingMatrix');
    
    AllStats{fileidx,1} = [{PeakTimeMatrix}, {PeakAmplitudeMatrix}, {PeakERA}, {TimeToPeakERA}, {FWHMMatrix}, {PeakOnsetMatrix}];
    SubjFittedHRFs{fileidx} = FittedHRFs;
    SubjFittedFullHRFs{fileidx} = FittedFullHRFs;
    SubjERAs{fileidx} = EventAvg;
    
    stattable(fileidx,1:3) = PeakTimeMatrix(:,2)';
    stattable(fileidx,4:6) = PeakAmplitudeMatrix(:,2)';
    stattable(fileidx, 7:9) = PeakOnsetMatrix(:,2)';
    stattable(fileidx,10) = PeakTimeMatrix(2,1);
    stattable(fileidx,11) = PeakAmplitudeMatrix(2,1);
    stattable(fileidx,12) = PeakOnsetMatrix(2,1);
    
end



%%
for subj = 1:numel(SubjFittedHRFs)
    for lengths=1:3
        for ints=2
            HRFmatrix{lengths,ints}(subj,:) = SubjFittedHRFs{subj}{lengths,ints};
            FullHRFmatrix{lengths,ints}(subj,:) = SubjFittedFullHRFs{subj}{lengths,ints};
            ERAmatrix{lengths,ints}(subj,:) = SubjERAs{subj}{lengths,ints}';
        end
    end
    for lengths=2
        for ints = 1
            HRFmatrix{lengths,ints}(subj,:) = SubjFittedHRFs{subj}{lengths,ints};
            FullHRFmatrix{lengths,ints}(subj,:) = SubjFittedFullHRFs{subj}{lengths,ints};
            ERAmatrix{lengths,ints}(subj,:) = SubjERAs{subj}{lengths,ints}';
        end
    end
end
%%
for lengths=1:3
    for ints=2
        HRFmatrix{lengths,ints} = smooth(nanmean(HRFmatrix{lengths, ints},1),0.1,'loess')';
        FullHRFmatrix{lengths,ints} = nanmean(FullHRFmatrix{lengths, ints},1);
        ERAmatrix{lengths,ints} = nanmean(ERAmatrix{lengths,ints},1);
    end
end

for lengths=2
    for ints=1
        HRFmatrix{lengths,ints} = smooth(nanmean(HRFmatrix{lengths, ints},1),0.1,'loess')';
        FullHRFmatrix{lengths,ints} = nanmean(FullHRFmatrix{lengths, ints},1);
        ERAmatrix{lengths,ints} = nanmean(ERAmatrix{lengths,ints},1);
    end
end

%% EFFECTS WITH WITHIN-SUBJ CORRECTIONS

TTP = table2array(stattable(:,1:3));
AMP = table2array(stattable(:,4:6));
ONS = table2array(stattable(:,7:9));

TTP = table2array(cat(2, stattable(:,1:3), stattable(:,10)));
AMP = table2array(cat(2, stattable(:,4:6), stattable(:,11)));
ONS = table2array(cat(2, stattable(:,7:9), stattable(:,12)));

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
%%
% Difference score analysis
TTPd = [TTP(:,3) - TTP(:,1), TTP(:,2) - TTP(:,1), TTP(:,3) - TTP(:,2), TTP(:,2) - TTP(:,4)];
AMPd = [AMP(:,3) - AMP(:,1), AMP(:,2) - AMP(:,1), AMP(:,3) - AMP(:,2), AMP(:,2) - AMP(:,4)];
ONSd = [ONS(:,3) - ONS(:,1), ONS(:,2) - ONS(:,1), ONS(:,3) - ONS(:,2), ONS(:,2) - ONS(:,4)];

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
TTPeffect(1,4) = TTPdstat(1,4) / sqrt( (nanvar(TTP(:,2),1) + nanvar(TTP(:,4),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));

AMPeffect(1,1) = AMPdstat(1,1) / sqrt( (nanvar(AMP(:,1),1) + nanvar(AMP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,2) = AMPdstat(1,2) / sqrt( (nanvar(AMP(:,1),1) + nanvar(AMP(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,3) = AMPdstat(1,3) / sqrt( (nanvar(AMP(:,2),1) + nanvar(AMP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,4) = AMPdstat(1,4) / sqrt( (nanvar(AMP(:,2),1) + nanvar(AMP(:,4),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));

ONSeffect(1,1) = ONSdstat(1,1) / sqrt( (nanvar(ONS(:,1),1) + nanvar(ONS(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,2) = ONSdstat(1,2) / sqrt( (nanvar(ONS(:,1),1) + nanvar(ONS(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,3) = ONSdstat(1,3) / sqrt( (nanvar(ONS(:,2),1) + nanvar(ONS(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,4) = ONSdstat(1,4) / sqrt( (nanvar(ONS(:,2),1) + nanvar(ONS(:,4),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));

TTPdstat(6,:) = TTPeffect;
AMPdstat(6,:) = AMPeffect;
ONSdstat(6,:) = ONSeffect;

writematrix(TTPstat, sprintf('%s', 'TimeToPeak', '_', 'durationManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'Amplitude', '_', 'durationManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'PeakOnset', '_', 'durationManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'TimeToPeakD', '_', 'durationManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'AmplitudeD', '_', 'durationManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'PeakOnsetD', '_', 'durationManipulation', '.csv'))

%% EFFECTS WITHOUT WITHIN-SUBJ CORRECTIONS

TTP = table2array(stattable(:,1:3));
AMP = table2array(stattable(:,4:6));
ONS = table2array(stattable(:,7:9));

TTP = table2array(cat(2, stattable(:,1:3), stattable(:,10)));
AMP = table2array(cat(2, stattable(:,4:6), stattable(:,11)));
ONS = table2array(cat(2, stattable(:,7:9), stattable(:,12)));

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
%%
% Difference score analysis
TTPd = [TTP(:,3) - TTP(:,1), TTP(:,2) - TTP(:,1), TTP(:,3) - TTP(:,2), TTP(:,2) - TTP(:,4)];
AMPd = [AMP(:,3) - AMP(:,1), AMP(:,2) - AMP(:,1), AMP(:,3) - AMP(:,2), AMP(:,2) - AMP(:,4)];
ONSd = [ONS(:,3) - ONS(:,1), ONS(:,2) - ONS(:,1), ONS(:,3) - ONS(:,2), ONS(:,2) - ONS(:,4)];

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
TTPeffect(1,4) = TTPdstat(1,4) / sqrt( (nanvar(TTP(:,2),1) + nanvar(TTP(:,4),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));

AMPeffect(1,1) = AMPdstat(1,1) / sqrt( (nanvar(AMP(:,1),1) + nanvar(AMP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,2) = AMPdstat(1,2) / sqrt( (nanvar(AMP(:,1),1) + nanvar(AMP(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,3) = AMPdstat(1,3) / sqrt( (nanvar(AMP(:,2),1) + nanvar(AMP(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
AMPeffect(1,4) = AMPdstat(1,4) / sqrt( (nanvar(AMP(:,2),1) + nanvar(AMP(:,4),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));

ONSeffect(1,1) = ONSdstat(1,1) / sqrt( (nanvar(ONS(:,1),1) + nanvar(ONS(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,2) = ONSdstat(1,2) / sqrt( (nanvar(ONS(:,1),1) + nanvar(ONS(:,2),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,3) = ONSdstat(1,3) / sqrt( (nanvar(ONS(:,2),1) + nanvar(ONS(:,3),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));
ONSeffect(1,4) = ONSdstat(1,4) / sqrt( (nanvar(ONS(:,2),1) + nanvar(ONS(:,4),1)) / 2) * ((size(TTPd,1)-3) / (size(TTPd,1)-2.25)) * (sqrt((size(TTPd,1)-2) / size(TTPd,1)));

TTPdstat(6,:) = TTPeffect;
AMPdstat(6,:) = AMPeffect;
ONSdstat(6,:) = ONSeffect;

writematrix(TTPstat, sprintf('%s', 'TimeToPeak', '_', 'durationManipulation', '.csv'))
writematrix(AMPstat, sprintf('%s', 'Amplitude', '_', 'durationManipulation', '.csv'))
writematrix(ONSstat, sprintf('%s', 'PeakOnset', '_', 'durationManipulation', '.csv'))
writematrix(TTPdstat, sprintf('%s', 'TimeToPeakD', '_', 'durationManipulation', '.csv'))
writematrix(AMPdstat, sprintf('%s', 'AmplitudeD', '_', 'durationManipulation', '.csv'))
writematrix(ONSdstat, sprintf('%s', 'PeakOnsetD', '_', 'durationManipulation', '.csv'))