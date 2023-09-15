clearvars();
cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis');
load('BigKahunaMatrix.mat');

subjectIDs = importdata('subject_list.txt');

for scount = 1:numel(subjectIDs)
    matchcount = 1;
    for rcount = 1:numel(BigKahunaMatrix)
        if BigKahunaMatrix(rcount).SID == subjectIDs{scount}(1:7)
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
        %if SecondSort{scount}{ucount}.ROI == 'M1.txt'
            %SecondSort{scount}{ucount} = [];
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
        CurrentV(datacount).VisualEvs = SubjectVisuals{scount}{datacount}.VSTARTSEC / 1 + 1;
        CurrentV(datacount).VisualEvs = round(CurrentV(datacount).VisualEvs, 0) - 4;
        CurrentV(datacount).VisualLengths = SubjectVisuals{scount}{datacount}.VLENGTH;
        CurrentV(datacount).VisualInts = SubjectVisuals{scount}{datacount}.VINTENS;
        for eventcount = 1:numel(CurrentV(datacount).VisualEvs)
        CurrentV(datacount).Timecourse{eventcount} = CurrentV(datacount).RawData(CurrentV(datacount).VisualEvs(eventcount)-5 : CurrentV(datacount).VisualEvs(eventcount)+19);
        end
        cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\PSubjData');
        save(sprintf('%s', ['V1_', subjectIDs{scount}(1:7)]), 'CurrentV');
    end
end

for scount = 1:numel(subjectIDs)
    for datacount = 1:numel(SubjectAuditory{scount})
        CurrentA(datacount).RawData = SubjectAuditory{scount}{datacount}.DATA;
        CurrentA(datacount).AuditoryEvs = SubjectAuditory{scount}{datacount}.ASTARTSEC / 1 + 1;
        CurrentA(datacount).AuditoryEvs = round(CurrentA(datacount).AuditoryEvs, 0) - 4;
        CurrentA(datacount).AuditoryLengths = SubjectAuditory{scount}{datacount}.ALENGTH;
        CurrentA(datacount).AuditoryInts = SubjectAuditory{scount}{datacount}.AINTENS;
        CurrentA(datacount).AuditoryInts = round(CurrentA(datacount).AuditoryInts,1);
        for eventcount = 1:numel(CurrentA(datacount).AuditoryEvs)-1
        CurrentA(datacount).Timecourse{eventcount} = CurrentA(datacount).RawData(CurrentA(datacount).AuditoryEvs(eventcount)-5 : CurrentA(datacount).AuditoryEvs(eventcount)+19);
        end
        cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\PSubjData');
        save(sprintf('%s', ['A1_', subjectIDs{scount}(1:7)]), 'CurrentA');
    end
end