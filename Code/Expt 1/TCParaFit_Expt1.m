clearvars();
cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\Expt1_GroupMatrix');
load('BigKahunaMatrix2.mat');

%% IV (condition) specification
eventtype = 'V';

conddef = {["Length 0.1 Int 0.1" "Length 0.1 Int 0.3" "Length 0.1 Int 0.9"; "Length 0.3 Int 0.1" "Length 0.3 Int 0.3" "Length 0.3 Int 0.9"; "Length 0.9 Int 0.1" "Length 0.9 Int 0.3" "Length 0.9 Int 0.9"];...
    ["Length 0.1 Int 0.01" "Length 0.1 Int 0.1" "Length 0.1 Int 1"; "Length 0.3 Int 0.01" "Length 0.3 Int 0.1" "Length 0.3 Int 1"; "Length 0.9 Int 0.01" "Length 0.9 Int 0.1" "Length 0.9 Int 1"]};

Amatcnt = 1;
Vmatcnt = 1;

for datacnt = 1:numel(BigKahunaMatrix)
    if BigKahunaMatrix(datacnt).ROI == "A1.txt"
        AData(Amatcnt) = BigKahunaMatrix(datacnt);
        Amatcnt = Amatcnt+1;
    elseif BigKahunaMatrix(datacnt).ROI == "V1.txt"
        VData(Vmatcnt) = BigKahunaMatrix(datacnt);
        Vmatcnt = Vmatcnt+1;
    end
end

%%
if eventtype == 'V'
    conddefinition = conddef{2};
    eventlengths = [0.1 0.3 0.9];
    eventintensities = [0.01 0.1 1];
    outputdir = "C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\out\vis_newscript"

    for datacnt = 1:numel(VData)
    VName(:,datacnt) = VData(datacnt).SID;
    VOnset(:,datacnt) = VData(datacnt).VSTARTSEC;
    VLength(:,datacnt) = VData(datacnt).VLENGTH;
    VInt(:,datacnt) = VData(datacnt).VINTENS;
    end
    VOnset = floor((VOnset / 0.625) - 4);

    for datacnt = 1:numel(VData)
    for evcnt = 1:size(VOnset(:,datacnt),1)
        VTC(:,evcnt,datacnt) = VData(datacnt).DATA((VOnset(evcnt,datacnt)-9):(VOnset(evcnt,datacnt)+31));
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

    ProtoMatrix{3, 3} = [];

    for Lcon = 1:numel(eventlengths)
    for Icon = 1:numel(eventintensities)
        ConditionMap{Lcon, Icon} = [eventlengths(Lcon) eventintensities(Icon)];
        for datacnt = 1:size(Lsearch{Lcon},2)
            matcnt = 1;
            for evcnt = 1:size(Lsearch{Lcon}(:,datacnt),1)
                if Lsearch{Lcon}(evcnt, datacnt) == 1 && Isearch{Icon}(evcnt, datacnt) == 1
                    ProtoMatrix{Lcon, Icon}(1:41, matcnt, datacnt) = VTC(1:end, evcnt,datacnt);
                    matcnt = matcnt+1;
                end
            end
        end
    end
    end

elseif eventtype == 'A'
    conddefinition = conddef{1};
    eventlengths = [0.1 0.3 0.9];
    eventintensities = [0.1 0.3 1];
    outputdir = "C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\out\aud_newscript"

    for datacnt = 1:numel(AData)
    AOnset(:,datacnt) = AData(datacnt).ASTARTSEC;
    ALength(:,datacnt) = AData(datacnt).ALENGTH;
    AInt(:,datacnt) = round(AData(datacnt).AINTENS, 1);
    end
    AOnset = floor((AOnset / 0.625) - 4);

    for datacnt = 1:numel(AData)
    for evcnt = 1:size(AOnset(:,datacnt),1)
        ATC(:,evcnt,datacnt) = AData(datacnt).DATA((AOnset(evcnt,datacnt)-9):(AOnset(evcnt,datacnt)+31));
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
    
    ProtoMatrix{3, 3} = [];

    for Lcon = 1:numel(eventlengths)
    for Icon = 1:numel(eventintensities)
        ConditionMap{Lcon, Icon} = [eventlengths(Lcon) eventintensities(Icon)];
        for datacnt = 1:size(Lsearch{Lcon},2)
            matcnt = 1;
            for evcnt = 1:size(Lsearch{Lcon}(:,datacnt),1)
                if Lsearch{Lcon}(evcnt, datacnt) == 1 && Isearch{Icon}(evcnt, datacnt) == 1
                    ProtoMatrix{Lcon, Icon}(1:41, matcnt, datacnt) = ATC(1:end, evcnt,datacnt);
                    matcnt = matcnt+1;
                end
            end
        end
    end
    end
    
end

for conC = 1:numel(ProtoMatrix)
    ProtoMatrix{conC} = reshape(ProtoMatrix{conC}, [41 size(ProtoMatrix{conC},2)*size(ProtoMatrix{conC},3)]);
end

for conC = 1:numel(ProtoMatrix)
    BaselineVector{conC} = ProtoMatrix{conC}(13, :);
    for i = 1:size(ProtoMatrix{conC}, 2)
        ProtoMatrix{conC}(:, i) = ((ProtoMatrix{conC}(:, i) - BaselineVector{conC}(:, i)) / BaselineVector{conC}(:, i)) * 100;
    end
    EventAvg{conC} = nanmean(ProtoMatrix{conC}(10:34,:), 2);
    PeakERA{conC} = max(EventAvg{conC});
    TimeToPeakERA{conC} = find(abs(EventAvg{conC}-PeakERA{conC}) < 0.001) * 0.625;
    ErrMargin{conC} = (nanstd(ProtoMatrix{conC}(10:34,:),0,2)/sqrt(size(ProtoMatrix{conC},2)))*1.96;
end

    EventAvg = reshape(EventAvg, [3 3]);
    ErrMargin = reshape(ErrMargin, [3 3]);

%% 
for lengths = 1:size(ProtoMatrix,1)
for ints = 1:size(ProtoMatrix,2)
    [FittingMatrix{lengths, ints}] = GammaFit(ProtoMatrix{lengths, ints}(10:34,:)', 0.625, 15);
    DelayMatrix{lengths, ints} = FittingMatrix{lengths, ints}.delay;
    AmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.amplitude;
    OnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.onset;
    DispersionMatrix{lengths, ints} = FittingMatrix{lengths, ints}.dispersion;
    PeakTimeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxIndexInSeconds;
    PeakAmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxValue;
    FWHMMatrix{lengths, ints} = FittingMatrix{lengths, ints}.FWHMInSeconds;
    PeakOnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.OnsetTimeInSeconds;
end
end

save('grouphrfdata.mat','FittingMatrix', 'PeakAmplitudeMatrix', 'PeakTimeMatrix', 'PeakOnsetMatrix');
save('grouperadata.mat','EventAvg', 'ErrMargin', 'PeakERA', 'TimeToPeakERA', 'ProtoMatrix');
