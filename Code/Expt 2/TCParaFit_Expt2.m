clearvars();
cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt2(2023)');
load('BigKahunaMatrix.mat');

addpath('C:\Users\alvinw\iCloudDrive\Documents\tron_analysis');
%% IV (condition) specification
eventtype = 'V';

conddef = {["NA" "Length 1.2 Int 0.3" "NA"; "Length 1.0 Int 0.9" "Length 1.2 Int 0.9" "Length 1.8 Int 0.9"]';...
    ["NA" "Length 1.2 Int 0.1" "NA"; "Length 1.0 Int 1" "Length 1.2 Int 1" "Length 1.8 Int 1"]'};

Amatcnt = 1;
Vmatcnt = 1;
%Mmatcnt = 1;

for datacnt = 1:numel(BigKahunaMatrix)
    if BigKahunaMatrix(datacnt).ROI == "A1.txt"
        AData(Amatcnt) = BigKahunaMatrix(datacnt);
        Amatcnt = Amatcnt+1;
    elseif BigKahunaMatrix(datacnt).ROI == "V1.txt"
        VData(Vmatcnt) = BigKahunaMatrix(datacnt);
        Vmatcnt = Vmatcnt+1;
%     elseif BigKahunaMatrix(datacnt).ROI == "M1"
%         MData(Mmatcnt) = BigKahunaMatrix(datacnt);
%         Mmatcnt = Mmatcnt+1;
    end
end

%%
if eventtype == 'V'
    conddefinition = conddef{2};
    eventlengths = [1.0 1.2 1.8];
    eventintensities = [0.1 1];
    
    for datacnt = 1:numel(VData)
    VName(:,datacnt) = VData(datacnt).SID;
    VOnset(:,datacnt) = VData(datacnt).VSTARTSEC;
    VLength(:,datacnt) = VData(datacnt).VLENGTH;
    VInt(:,datacnt) = VData(datacnt).VINTENS;
    end
    VOnset = floor((VOnset / 1) - 4);

    for datacnt = 1:numel(VData)
    for evcnt = 1:size(VOnset(:,datacnt),1)
        VTC(:,evcnt,datacnt) = VData(datacnt).DATA((VOnset(evcnt,datacnt)-4):(VOnset(evcnt,datacnt)+15));
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

elseif eventtype == 'A'
    conddefinition = conddef{1};
    eventlengths = [1.0 1.2 1.8];
    eventintensities = [0.3 1];

    for datacnt = 1:numel(AData)
    AOnset(:,datacnt) = AData(datacnt).ASTARTSEC;
    ALength(:,datacnt) = AData(datacnt).ALENGTH;
    AInt(:,datacnt) = round(AData(datacnt).AINTENS, 1);
    end
    AOnset = floor((AOnset / 1) - 4);

    for datacnt = 1:numel(AData)
    for evcnt = 1:size(AOnset(:,datacnt),1)
        ATC(:,evcnt,datacnt) = AData(datacnt).DATA((AOnset(evcnt,datacnt)-4):(AOnset(evcnt,datacnt)+15));
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

for conC = 1:numel(ProtoMatrix)
    ProtoMatrix{conC} = reshape(ProtoMatrix{conC}, [20, size(ProtoMatrix{conC},2)*size(ProtoMatrix{conC},3)]);
end

for conC = 1:numel(ProtoMatrix)
    BaselineVector{conC} = ProtoMatrix{conC}(4, :);
    for i = 1:size(ProtoMatrix{conC}, 2)
        ProtoMatrix{conC}(:, i) = ((ProtoMatrix{conC}(:, i) - BaselineVector{conC}(:, i)) / BaselineVector{conC}(:, i)) * 100;
    end
    EventAvg{conC} = nanmean(ProtoMatrix{conC}(5:20,:), 2);
    PeakERA{conC} = max(EventAvg{conC});
    TimeToPeakERA{conC} = find(abs(EventAvg{conC}-PeakERA{conC}) < 0.001) * 1;
    ErrMargin{conC} = (nanstd(ProtoMatrix{conC}(5:20,:),0,2)/sqrt(size(ProtoMatrix{conC},2)))*1.96;
end

    EventAvg = reshape(EventAvg, [3 2]);
    ErrMargin = reshape(ErrMargin, [3 2]);
%     
%     cd('C:\Users\Alvin\iCloudDrive\TRONAnalysis2020')
% 
%     save(sprintf('%s', [eventtype, '_ERAresults']), 'ProtoMatrix', 'EventAvg', 'ErrMargin');

% for condition = 1:numel(XtremeProtoMatrix)
%     plotTitle = sprintf('%s', [conddefinition(condition)]);
%     TCPlot1(looks, trvector, EventAvg{condition}, ErrMargin{condition}, plotTitle, conddefinition(condition));
% end

%% 
for lengths = 1:size(ProtoMatrix,1)
for ints = 2 %1:size(ProtoMatrix,2)

    [FittingMatrix{lengths, ints}] = GammaFit(ProtoMatrix{lengths, ints}(5:20,:)', 1, 15);
    DelayMatrix{lengths, ints} = FittingMatrix{lengths, ints}.delay;
    AmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.amplitude;
    OnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.onset;
    DispersionMatrix{lengths, ints} = FittingMatrix{lengths, ints}.dispersion;
    PeakTimeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxIndexInSeconds;
    PeakAmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxValue;
    GammaParams = [FittingMatrix{lengths, ints}.amplitude FittingMatrix{lengths, ints}.delay FittingMatrix{lengths, ints}.onset FittingMatrix{lengths, ints}.dispersion FittingMatrix{lengths, ints}.MaxIndexInSeconds FittingMatrix{lengths, ints}.MaxValue];
    
    save(sprintf('%s', [eventtype, '_', 'FittingMatrix', '.mat']), 'FittingMatrix');
    save(sprintf('%s', [eventtype, '_', '_GammaParams_', 'l', num2str(lengths), '_i', num2str(ints), '.mat']),'GammaParams');
    hold on;
    plot(0:1:15, nanmean(ProtoMatrix{lengths, ints}(5:20, :), 2), 'r', 'LineWidth', 2)
    
    savefig(sprintf('%s', [eventtype, '_SingleGammaFit_', 'l', num2str(lengths), '_i', num2str(ints), '.fig']));
end
end
    
hold off;

%%
linecolor = ['r' 'g' 'b'];
k = 1;
for ints = 2
    for lengths = 1:3
        timecoursePlot(lengths) = plot(0:1:15, nanmean(ProtoMatrix{lengths, ints}(5:20, :),2), linecolor(lengths), 'LineWidth', 1, 'DisplayName', sprintf('%s', conddefinition(lengths, ints)));
        hold on;
        TR = 1/1000; % Pretend TR is 1/1000 of actual for upsampling purposes
        eventlength = 15*1000;
        t = 0:TR:15-TR;

        plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(lengths), 'LineWidth', 2);
        hold on;
        plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(lengths), 'LineWidth', 2);
        plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(lengths), 'LineWidth', 2);
        hold on;
    end
    
    plotTitle = sprintf('%s', ["Comparison of Gamma Function Fits for ", eventtype, "1 ", "by stimulus intensity for duration ", eventlengths(lengths)]);
    xlabel('Time from event onset (secs)', 'fontsize', 14); % x-axis label
    ylabel('Response', 'fontsize', 14); % y-axis label
    ylim([-1 2]);
    set(gca, 'fontsize', 14);
    title([plotTitle],'fontsize',14); % plot title

    legend(timecoursePlot(1:3), [sprintf('%s', conddefinition(1, 2)), '; a: ' num2str(PeakAmplitudeMatrix{1, ints}) ', d: ' num2str(PeakTimeMatrix{1, ints})],...
    [sprintf('%s', conddefinition(2, 2)), '; a: ' num2str(PeakAmplitudeMatrix{2, ints}) ', d: ' num2str(PeakTimeMatrix{2, ints})],...
    [sprintf('%s', conddefinition(3, 2)), '; a: ' num2str(PeakAmplitudeMatrix{3, ints}) ', d: ' num2str(PeakTimeMatrix{3, ints})]);

    %cd('/Volumes/Passport/TRoN/Analysis/PGroup2');
    savefig(sprintf('%s', [eventtype, '_GammaFitLengthsV2_', 'l', num2str(k), '_i', num2str(ints), '.fig']));
    hold off;
    k = k+1;
end
%%
for lengths = 2%1:size(ProtoMatrix,1)
for ints = 1%:size(ProtoMatrix,2)

    [FittingMatrix{lengths, ints}] = GammaFit(ProtoMatrix{lengths, ints}(5:20,:)', 1, 15);
    DelayMatrix{lengths, ints} = FittingMatrix{lengths, ints}.delay;
    AmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.amplitude;
    OnsetMatrix{lengths, ints} = FittingMatrix{lengths, ints}.onset;
    DispersionMatrix{lengths, ints} = FittingMatrix{lengths, ints}.dispersion;
    PeakTimeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxIndexInSeconds;
    PeakAmplitudeMatrix{lengths, ints} = FittingMatrix{lengths, ints}.MaxValue;
    GammaParams = [FittingMatrix{lengths, ints}.amplitude FittingMatrix{lengths, ints}.delay FittingMatrix{lengths, ints}.onset FittingMatrix{lengths, ints}.dispersion FittingMatrix{lengths, ints}.MaxIndexInSeconds FittingMatrix{lengths, ints}.MaxValue];
    
    save(sprintf('%s', [eventtype, '_', 'FittingMatrixI', '.mat']), 'FittingMatrix');
    save(sprintf('%s', [eventtype, '_', '_GammaParams_', 'l', num2str(lengths), '_i', num2str(ints), '.mat']),'GammaParams');
    hold on;
    plot(0:1:15, nanmean(ProtoMatrix{lengths, ints}(5:20, :), 2), 'r', 'LineWidth', 2)
    
    savefig(sprintf('%s', [eventtype, '_SingleGammaFit_', 'l', num2str(lengths), '_i', num2str(ints), '.fig']));
end
end
    
hold off;

%%
linecolor = ['r' 'g' 'b'];
k = 1;
for ints = 1:2
    for lengths = 2
        timecoursePlot2(ints) = plot(0:1:15, nanmean(ProtoMatrix{lengths, ints}(5:20, :),2), linecolor(ints), 'LineWidth', 1, 'DisplayName', sprintf('%s', conddefinition(lengths, ints)));
        hold on;
        TR = 1/1000; % Pretend TR is 1/1000 of actual for upsampling purposes
        eventlength = 15*1000;
        t = 0:TR:15-TR;

        plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(ints), 'LineWidth', 2);
        hold on;
        plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(ints), 'LineWidth', 2);
        %plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(lengths), 'LineWidth', 2);
        hold on;
    end
    
    plotTitle = sprintf('%s', ["Comparison of Gamma Function Fits for ", eventtype, "1 ", "by stimulus intensity for duration ", eventlengths(lengths)]);
    xlabel('Time from event onset (secs)', 'fontsize', 14); % x-axis label
    ylabel('Response', 'fontsize', 14); % y-axis label
    ylim([-1 2]);
    set(gca, 'fontsize', 14);
    title([plotTitle],'fontsize',14); % plot title

    
%     legend(timecoursePlot(1:3), [sprintf('%s', conddefinition(1, 1)), '; a: ' num2str(PeakAmplitudeMatrix{1, ints}) ', d: ' num2str(PeakTimeMatrix{1, ints})],...
%     [sprintf('%s', conddefinition(2, 1)), '; a: ' num2str(PeakAmplitudeMatrix{2, ints}) ', d: ' num2str(PeakTimeMatrix{2, ints})],...
%     [sprintf('%s', conddefinition(3, 1)), '; a: ' num2str(PeakAmplitudeMatrix{3, ints}) ', d: ' num2str(PeakTimeMatrix{3, ints})]);

    %cd('/Volumes/Passport/TRoN/Analysis/PGroup2');
    
end
legend(timecoursePlot2(1:2), [sprintf('%s', conddefinition(2, 1)), '; a: ' num2str(PeakAmplitudeMatrix{2, 1}) ', d: ' num2str(PeakTimeMatrix{2, 1})],...
    [sprintf('%s', conddefinition(2, 2)), '; a: ' num2str(PeakAmplitudeMatrix{2, 2}) ', d: ' num2str(PeakTimeMatrix{2, 2})]);
savefig(sprintf('%s', [eventtype, '_GammaFitLengthsV2_', 'i', num2str(k), '_l', num2str(ints), '.fig']));
    hold off;
    k = k+1;
%%
k=1;
% for ints = 1:3
%     for lengths = 1:size(ProtoMatrix,1)
% %         results.amplitude = nanmean(AmplitudeMatrix{lengths, ints});
% %         results.delay = nanmean(DelayMatrix{lengths, ints});
% %         results.onset = nanmean(OnsetMatrix{lengths, ints});
% %         results.dispersion = nanmean(DispersionMatrix{lengths, ints});
% %         MoreResults{lengths, ints} = CanonicalGammaInfo(results, 1, 15);

%         timecoursePlot(lengths, ints) = plot(0:1:15, nanmean(ProtoMatrix{lengths, ints}(11:35, :),2), linecolor(lengths), 'LineWidth', 1, 'DisplayName', sprintf('%s', conddefinition(lengths, ints)));
%         hold on;

%         TR = 1/1000; % Pretend TR is 1/1000 of actual for upsampling purposes
%         eventlength = 15*1000;
%         t = 0:TR:15-TR;
%         plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(lengths), 'LineWidth', 2);
%         hold on;
%         plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(lengths), 'LineWidth', 2);
%         plot(t,FittingMatrix{lengths, ints}.hrf(1,:), linecolor(lengths), 'LineWidth', 2);
%         hold on;
%     end
    
%     plotTitle = sprintf('%s', ["Comparison of Gamma Function Fits for ", eventtype, "1 ", "by stimulus duration for intensity ", eventlengths(ints)]);
%     xlabel('Time from event onset (secs)', 'fontsize', 14); % x-axis label
%     ylabel('Response', 'fontsize', 14); % y-axis label
%     ylim([-1 2]);
%     set(gca, 'fontsize', 14);
%     title([plotTitle],'fontsize',14); % plot title

%     legend(timecoursePlot(1:3, ints), [sprintf('%s', conddefinition(1, ints)), '; a: ' num2str(PeakAmplitudeMatrix{1, ints}) ', d: ' num2str(PeakTimeMatrix{1, ints})],...
%     [sprintf('%s', conddefinition(2, ints)), '; a: ' num2str(PeakAmplitudeMatrix{2, ints}) ', d: ' num2str(PeakTimeMatrix{2, ints})],...
%     [sprintf('%s', conddefinition(3, ints)), '; a: ' num2str(PeakAmplitudeMatrix{3, ints}) ', d: ' num2str(PeakTimeMatrix{3, ints})]);
%     %cd('/Volumes/Passport/TRoN/Analysis/PGroup2');
%     savefig(sprintf('%s', [eventtype, '_GammaFitIntsV2_', 'l', num2str(k), '_i', num2str(ints), '.fig']));
%     hold off;
%     k = k+1;
% end
