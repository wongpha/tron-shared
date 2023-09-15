clearvars();

% Per-subject ROI timecourse event-related averaging, linear gamma function
% fitting and point estimate computation

% Written Sep 2018 by Alvin Wong for TRoN Project

clearvars();
cd('/Volumes/Passport/TRoN/TRoN_data');
load('BigKahunaMatrix.mat');

subjectIDs = importdata('subjectIDs.txt');

for scount = 1:numel(subjectIDs)
    matchcount = 1;
    for rcount = 1:numel(BigKahunaMatrix)
        if BigKahunaMatrix(rcount).SID == subjectIDs{scount}
            FirstSort{matchcount, scount} = BigKahunaMatrix(rcount);
            matchcount = matchcount + 1;
        end
    end
end

for scount = 1:numel(subjectIDs)
    matchcountV = 1;
    matchcountA = 1;
    for tcount = 1:numel(FirstSort(:,scount))
        if isempty(FirstSort{tcount,scount}) ~= 1
            SecondSort{scount}{tcount} = FirstSort{tcount,scount};
        end
    end
    for ucount = 1:numel(SecondSort{scount})
        if SecondSort{scount}{ucount}.ROI == "V1.txt"
            SubjectVisuals{scount}{matchcountV} = SecondSort{scount}{ucount};
            matchcountV = matchcountV+1;
        elseif SecondSort{scount}{ucount}.ROI == "A1.txt"
            SubjectAuditory{scount}{matchcountA} = SecondSort{scount}{ucount};
            matchcountA = matchcountA+1;
        end
    end
end
%%
for scount = 1:numel(subjectIDs)
    for datacount = 1:numel(SubjectVisuals{scount})
        CurrentV(datacount).RawData = SubjectVisuals{scount}{datacount}.DATA;
        CurrentV(datacount).VisualEvs = SubjectVisuals{scount}{datacount}.VSTARTSEC / 0.625 + 1;
        CurrentV(datacount).VisualEvs = round(CurrentV(datacount).VisualEvs, 0) - 4;
        CurrentV(datacount).VisualLengths = SubjectVisuals{scount}{datacount}.VLENGTH;
        CurrentV(datacount).VisualInts = SubjectVisuals{scount}{datacount}.VINTENS;
        cd('/Volumes/Passport/TRoN/Analysis/PSubj2/Data');
        save(sprintf('%s', ['V1_', subjectIDs{scount}]), 'CurrentV');
    end
end

for scount = 1:numel(subjectIDs)
    for datacount = 1:numel(SubjectAuditory{scount})
        CurrentA(datacount).RawData = SubjectAuditory{scount}{datacount}.DATA;
        CurrentA(datacount).AuditoryEvs = SubjectAuditory{scount}{datacount}.ASTARTSEC / 0.625 + 1;
        CurrentA(datacount).AuditoryEvs = round(CurrentA(datacount).AuditoryEvs, 0) - 4;
        CurrentA(datacount).AuditoryLengths = SubjectAuditory{scount}{datacount}.ALENGTH;
        CurrentA(datacount).AuditoryInts = SubjectAuditory{scount}{datacount}.AINTENS;
        CurrentA(datacount).AuditoryInts = round(CurrentA(datacount).AuditoryInts,1);
        cd('/Volumes/Passport/TRoN/Analysis/PSubj2/Data');
        save(sprintf('%s', ['A1_', subjectIDs{scount}]), 'CurrentA');
    end
end

cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\Data'); % Path to data matrices output by PerSubjProc
addpath 'C:\Users\Alvin\iCloudDrive\tron_analysis'
eventtype = "A"; % Set ROI type (A or V)
filedir = dir('A1*.mat'); % Set filename search criteria

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
                
    if eventtype == 'V' % V1 timecourse processing to build protomatrix
        conddefinition = conddef{2};
        eventlengths = [0.1 0.3 0.9];
        eventintensities = [0.01 0.1 1];
        roiname = 'V1';
        
        for datacnt = 1:numel(CurrentV)
        VOnset(:, datacnt) = CurrentV(datacnt).VisualEvs(1:end);
        VLength(:,datacnt) = CurrentV(datacnt).VisualLengths;
        VInt(:,datacnt) = CurrentV(datacnt).VisualInts;
        end
        
        for datacnt = 1:numel(CurrentV)
        for evcnt = 1:size(VOnset(:,datacnt),1)
            VTC(:,evcnt,datacnt) = CurrentV(datacnt).RawData((VOnset(evcnt,datacnt)-9):(VOnset(evcnt,datacnt)+31));
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

    elseif eventtype == 'A' % A1 timecourse processing to build protomatrix
        conddefinition = conddef{1};
        eventlengths = [0.1 0.3 0.9];
        eventintensities = [0.1 0.3 1];
        roiname = 'A1';
        
        for datacnt = 1:numel(CurrentA)
        AOnset(:, datacnt) = CurrentA(datacnt).AuditoryEvs(1:end);
        ALength(:,datacnt) = CurrentA(datacnt).AuditoryLengths;
        AInt(:,datacnt) = CurrentA(datacnt).AuditoryInts;
        end
        
        for datacnt = 1:numel(CurrentA)
        for evcnt = 1:size(AOnset(:,datacnt),1)
            if AOnset(evcnt,datacnt) == 590
                AOnset(evcnt,datacnt) = AOnset(evcnt,datacnt) - 1;
            end
            ATC(:,evcnt,datacnt) = CurrentA(datacnt).RawData((AOnset(evcnt,datacnt)-9):(AOnset(evcnt,datacnt)+31));
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
    
    for conC = 1:numel(ProtoMatrix) % convert 3D cells in protomatrix to 2D cells
        ProtoMatrix{conC} = reshape(ProtoMatrix{conC}, [41 size(ProtoMatrix{conC},2)*size(ProtoMatrix{conC},3)]);
    end
    if eventtype == 'V'
        save(sprintf('%s', [roiname, '_', subjname, '_data']), 'ProtoMatrix', 'ConditionMap', 'VInt', 'VLength', 'VOnset', 'VTC');
    elseif eventtype == 'A'
        save(sprintf('%s', [roiname, '_', subjname, '_data']), 'ProtoMatrix', 'ConditionMap', 'AInt', 'ALength', 'AOnset', 'ATC');
    end

end
