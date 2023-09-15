% Figure builder
clearvars();
cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\out\aud_newscript');
erafiledir = dir('*ERAvgs*.mat');
fitfiledir = dir('*FittingMatrix.mat');
cd(pwd)
lengths = [100 300 900];
intensities = [10 30 90];
xticklabs = {'-1', '1', '3', '5', '7', '9', '11', '13', '15'};
        
for erafileidx = 1:numel(erafiledir)
    erafilename = sprintf('%s', [erafiledir(erafileidx).folder, filesep, erafiledir(erafileidx).name]);
    eradata{erafileidx,1} = load(erafilename);
    fitfilename = sprintf('%s', [fitfiledir(erafileidx).folder, filesep, fitfiledir(erafileidx).name]);
    fitdata{erafileidx,1} = load(fitfilename);
end

grouphrfdata = load('LoadThis.mat');
grouperadata = load('LoadThis2.mat');

%% OVERALL HRFs

for l = 1:numel(lengths)
    for i = 1:numel(intensities)
        for s = 1:numel(erafiledir)
        %eravecs{l,i}(s,:) = eradata{s}.EventAvg{l,i}
        nbpmagdist{l,i}(s,:) = fitdata{s}.FittingMatrix{l,i}.MaxValue;
        nbpdeldist{l,i}(s,:) = fitdata{s}.FittingMatrix{l,i}.MaxIndexInSeconds;
        end
        
        hrfvecs{l,i}= grouphrfdata.FittingMatrix{l,i}.hrf;
        hrfvec2{l,i} = grouphrfdata.FittingMatrix{l,i}.fullhrf;
        eravecs{l,i} = grouperadata.EventAvg{l,i};
        smoothhrfvec2{l,i} = hrfvec2{l,i};
        samplingRateIncrease = 1000;
        newSamplePoints = linspace(0,15,15000);
        smoothhrfvec2{l,i} = spline(0:0.625:15, smoothhrfvec2{l,i}, newSamplePoints);
    end
end


ybounds = [-0.2 0.8];
subplot(2,3,4)
for l = 1
    for i = 1:numel(intensities)
        cols = ['r' 'g' 'b'];
        plot(0:0.625:15, eravecs{l,i}, cols(i))
        hold on
        TR = 0.625/1000;
        eventlength = 15;
        t = 0:TR:eventlength-TR;
        t2 = 0:0.001:15-0.001;
        plot(t2, smoothhrfvec2{l,i}, cols(i),'LineWidth', 2)
        %plot(t, hrfvecs{l,i}, cols(i), 'LineWidth', 2)
        ylim(ybounds);
        yticks(ybounds(1):0.2:ybounds(2));
        xticks(0:2:16);
        xticklabels(xticklabs);
        ax = gca;
        ax.FontSize = 12; 
        ylabel('% signal change', 'fontsize', 15);
        title(sprintf('%s', [num2str(lengths(l)), 'ms']),'fontsize',16);
    end
end

subplot(2,3,5)
for l = 2
    for i = 1:numel(intensities)
        cols = ['r' 'g' 'b'];
        plot(0:0.625:15, eravecs{l,i}, cols(i))
        hold on
        TR = 0.625/1000;
        eventlength = 15;
        t = 0:TR:eventlength-TR;
        t2 = 0:0.001:15-0.001;
        plot(t2, smoothhrfvec2{l,i}, cols(i),'LineWidth', 2)
        %plot(t,hrfvecs{l,i}, cols(i), 'LineWidth', 2)
        ylim(ybounds);
        yticks(ybounds(1):0.2:ybounds(2));
        xticks(0:2:16);
        xticklabels(xticklabs);
        ax = gca;
        ax.FontSize = 12; 
        title(sprintf('%s', [num2str(lengths(l)), 'ms']),'fontsize',16);
    end
end

subplot(2,3,6)
for l = 3
    for i = 1:numel(intensities)
        cols = ['r' 'g' 'b'];
        plot(0:0.625:15, eravecs{l,i}, cols(i),'HandleVisibility','off')
        hold on
        TR = 0.625/1000;
        eventlength = 15;
        t = 0:TR:eventlength-TR;
        t2 = 0:0.001:15-0.001;
        plot(t2, smoothhrfvec2{l,i}, cols(i),'LineWidth', 2)
        %plot(t,hrfvecs{l,i}, cols(i), 'LineWidth', 2)
        ylim(ybounds);
        yticks(ybounds(1):0.2:ybounds(2));
        xticks(0:2:16);
        xticklabels(xticklabs);
        ax = gca;
        ax.FontSize = 12; 
        title(sprintf('%s', [num2str(lengths(l)), 'ms']),'fontsize',16);
        legend('10%', '30%', '90%','fontsize', 15);
    end
end

subplot(2,3,1)
for l = 1:numel(lengths)
    for i = 1
        cols = {[0.9 0.3 0.6],[0.4 0.9 0.9],[0.5 0.4 1]};
        plot(0:0.625:15, eravecs{l,i}, 'Color', cols{l})
        hold on
        TR = 0.625/1000;
        eventlength = 15;
        t = 0:TR:eventlength-TR;
        t2 = 0:0.001:15-0.001;
        plot(t2, smoothhrfvec2{l,i}, 'Color', cols{l},'LineWidth', 2)
        %plot(t, hrfvecs{l,i}, 'Color', cols{l}, 'LineWidth', 2)
        ylim(ybounds);
        
        yticks(ybounds(1):0.2:ybounds(2));
        xticks(0:2:16);
        xticklabels(xticklabs);
        ax = gca;
        ax.FontSize = 12; 
        ylabel('% signal change','fontsize', 15);
        title(sprintf('%s', [num2str(intensities(i)), '% intensity']),'fontsize',16);
    end
end

subplot(2,3,2)
for l = 1:numel(lengths)
    for i = 2
        cols = {[0.9 0.3 0.6],[0.4 0.9 0.9],[0.5 0.4 1]};
        plot(0:0.625:15, eravecs{l,i}, 'Color', cols{l})
        hold on
        TR = 0.625/1000;
        eventlength = 15;
        t = 0:TR:eventlength-TR;
        t2 = 0:0.001:15-0.001;
        plot(t2, smoothhrfvec2{l,i}, 'Color', cols{l},'LineWidth', 2)
        %plot(t,hrfvecs{l,i}, 'Color', cols{l}, 'LineWidth', 2)
        ylim(ybounds);
        yticks(ybounds(1):0.2:ybounds(2));
        xticks(0:2:16);
        xticklabels(xticklabs);
        ax = gca;
        ax.FontSize = 12; 
        xlabel('Time in seconds', 'fontsize', 15);
        title(sprintf('%s', [num2str(intensities(i)), '% intensity']),'fontsize',16);
    end
end

subplot(2,3,3)
for l = 1:numel(lengths)
    for i = 3
        cols = {[0.9 0.3 0.6],[0.4 0.9 0.9],[0.5 0.4 1]};
        plot(0:0.625:15, eravecs{l,i}, 'Color', cols{l}, 'HandleVisibility','off')
        hold on
        TR = 0.625/1000;
        eventlength = 15;
        t = 0:TR:eventlength-TR;
        t2 = 0:0.001:15-0.001;
        plot(t2, smoothhrfvec2{l,i}, 'Color', cols{l},'LineWidth', 2)
        %plot(t,hrfvecs{l,i},'Color', cols{l}, 'LineWidth', 2)
        ylim(ybounds);
        yticks(ybounds(1):0.2:ybounds(2));
        xticks(0:2:16);
        xticklabels(xticklabs);
        ax = gca;
        ax.FontSize = 12; 
        title(sprintf('%s', [num2str(intensities(i)), '% intensity']),'fontsize',16);
        legend('100ms', '300ms', '900ms', 'fontsize', 15);
    end
end


%% PEAK MAGNITUDE
nbstats = csvread('Amplitude_durationManipulation.csv',0,0);
bstats = csvread("C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out\auditory\BootStp_Amplitude_durationManipulation.csv",0,0);

figure;
subplot(1,3,1)

xaxes = [100 300 900];
means = plot(xaxes, nbstats(1,1:3), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,0+i) nbstats(5,0+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(10:12,3), 'kx');
%ci2 = plot(xaxes, nbstats(10:12,4), 'kx');

means = plot(xaxes, bstats(1,1:3), '-o', 'Color', [0.9 0.3 0.6], 'LineWidth', 2, 'MarkerFaceColor', [0.9 0.3 0.6], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)],  [bstats(5,0+i) bstats(6,0+i)],':', 'LineWidth', 3, 'Color', [0.9 0.3 0.6]);
end
%ci1 = plot(xaxes, bstats(10:12,3), 'o', 'Color', [0.9 0.3 0.6], 'MarkerFaceColor', [0.9 0.3 0.6]);
%ci2 = plot(xaxes, bstats(10:12,4), 'o', 'Color', [0.9 0.3 0.6], 'MarkerFaceColor', [0.9 0.3 0.6]);

xlim([0 1000]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
yticks([0:0.1:0.8]);
ylim([0 0.8]);
ax = gca;
ax.FontSize = 13; 
ylabel('Peak magnitude (% signal change)', 'fontsize', 15);

title(sprintf('%s', [num2str(10), '% intensity']),'fontsize',16);

subplot(1,3,2)
means = plot(xaxes, nbstats(1,4:6), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,3+i) nbstats(5,3+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(13:15,3), 'k.');
%ci2 = plot(xaxes, nbstats(13:15,4), 'k.');

means = plot(xaxes, bstats(1,4:6), '-o', 'Color', [0.4 0.9 0.9], 'LineWidth', 2, 'MarkerFaceColor', [0.4 0.9 0.9],  'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,3+i) bstats(6,3+i)], ':', 'LineWidth', 3, 'Color', [0.4 0.9 0.9]);
end
%ci1 = plot(xaxes, bstats(13:15,3), 'o', 'Color', [0.4 0.9 0.9],'MarkerFaceColor', [0.4 0.9 0.9]);
%ci2 = plot(xaxes, bstats(13:15,4), 'o', 'Color', [0.4 0.9 0.9],'MarkerFaceColor', [0.4 0.9 0.9]);

xlim([0 1000]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
yticks([0:0.1:0.8]);
ylim([0 0.8]);

ax = gca;
ax.FontSize = 13; 
xlabel('Stimulus duration (ms)', 'fontsize', 15);
title(sprintf('%s', [num2str(30), '% intensity']),'fontsize',16);

subplot(1,3,3)
means = plot(xaxes, nbstats(1,7:9), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,6+i) nbstats(5,6+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(16:18,3), 'k.');
%ci2 = plot(xaxes, nbstats(16:18,4), 'k.');

means = plot(xaxes, bstats(1,7:9), '-o', 'Color', [0.5 0.4 1], 'LineWidth', 2, 'MarkerFaceColor', [0.5 0.4 1], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,6+i) bstats(6,6+i)], ':', 'LineWidth', 3, 'Color', [0.5 0.4 1]);
end
%ci1 = plot(xaxes, bstats(16:18,3), 'o', 'Color', [0.5 0.4 1], 'MarkerFaceColor', [0.5 0.4 1]);
%ci2 = plot(xaxes, bstats(16:18,4), 'o', 'Color', [0.5 0.4 1], 'MarkerFaceColor', [0.5 0.4 1]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
xlim([0 1000]);
yticks([0:0.1:0.8]);
ylim([0 0.8]);
ax = gca;
ax.FontSize = 13; 
title(sprintf('%s', [num2str(90), '% intensity']),'fontsize',16);

hold off

%% TIME TO PEAK
nbstats = csvread('TimeToPeak_durationManipulation.csv',0,0);
bstats = csvread("C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out\auditory\BootStp_TimeToPeak_durationManipulation.csv",0,0);

figure;
subplot(2,3,1)

xaxes = [100 300 900];
means = plot(xaxes, nbstats(1,1:3), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,0+i) nbstats(5,0+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(10:12,7), 'k.');
%ci2 = plot(xaxes, nbstats(10:12,8), 'k.');

means = plot(xaxes, bstats(1,1:3), '-o', 'Color', [0.9 0.3 0.6], 'LineWidth', 2, 'MarkerFaceColor', [0.9 0.3 0.6], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,0+i) bstats(6,0+i)], ':', 'LineWidth', 3, 'Color', [0.9 0.3 0.6]);
end
%ci1 = plot(xaxes, bstats(10:12,7), 'o', 'Color', [0.9 0.3 0.6], 'MarkerFaceColor', [0.9 0.3 0.6]);
%ci2 = plot(xaxes, bstats(10:12,8), 'o', 'Color', [0.9 0.3 0.6], 'MarkerFaceColor', [0.9 0.3 0.6]);

xlim([0 1000]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
yticks([2.8:0.4:5.6]);
ylim([2.8 5.6]);
ax = gca;
ax.FontSize = 13; 
ylabel('Peak latency (s)', 'fontsize', 15);

title(sprintf('%s', [num2str(10), '% intensity']),'fontsize',16);

subplot(2,3,2)

xaxes = [100 300 900];
means = plot(xaxes, nbstats(1,4:6), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,3+i) nbstats(5,3+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(13:15,7), 'k.');
%ci2 = plot(xaxes, nbstats(13:15,8), 'k.');

means = plot(xaxes, bstats(1,4:6), '-o', 'Color', [0.4 0.9 0.9], 'LineWidth', 2, 'MarkerFaceColor', [0.4 0.9 0.9], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,3+i) bstats(6,3+i)],  ':','LineWidth', 3, 'Color', [0.4 0.9 0.9]);
end
%ci1 = plot(xaxes, bstats(13:15,7), 'o', 'Color', [0.4 0.9 0.9], 'MarkerFaceColor', [0.4 0.9 0.9]);
%ci2 = plot(xaxes, bstats(13:15,8), 'o', 'Color', [0.4 0.9 0.9], 'MarkerFaceColor', [0.4 0.9 0.9]);

xlim([0 1000]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
yticks([2.8:0.4:5.6]);
ylim([2.8 5.6]);
ax = gca;
ax.FontSize = 13; 
xlabel('Stimulus duration (ms)', 'fontsize', 15);
title(sprintf('%s', [num2str(30), '% intensity']),'fontsize',16);

subplot(2,3,3)

means = plot(xaxes, nbstats(1,7:9), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,6+i) nbstats(5,6+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(16:18,7), 'k.');
%ci2 = plot(xaxes, nbstats(16:18,8), 'k.');

means = plot(xaxes, bstats(1,7:9), '-o', 'Color', [0.5 0.4 1], 'LineWidth', 2, 'MarkerFaceColor', [0.5 0.4 1], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,6+i) bstats(6,6+i)],  ':','LineWidth', 3, 'Color', [0.5 0.4 1]);
end
%ci1 = plot(xaxes, bstats(16:18,7), 'o', 'Color', [0.5 0.4 1],'MarkerFaceColor', [0.5 0.4 1]);
%ci2 = plot(xaxes, bstats(16:18,8), 'o', 'Color', [0.5 0.4 1], 'MarkerFaceColor', [0.5 0.4 1]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
xlim([0 1000]);
yticks([2.8:0.4:5.6]);
ylim([2.8 5.6]);
ax = gca;
ax.FontSize = 13; 
title(sprintf('%s', [num2str(90), '% intensity']),'fontsize',16);

nbstats = csvread('PeakOnset_durationManipulation.csv',0,0);
bstats = csvread("C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out\auditory\BootStp_PeakOnset_durationManipulation.csv",0,0);

subplot(2,3,4)

xaxes = [100 300 900];
means = plot(xaxes, nbstats(1,1:3), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,0+i) nbstats(5,0+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(10:12,23), 'k.');
%ci2 = plot(xaxes, nbstats(10:12,24), 'k.');

means = plot(xaxes, bstats(1,1:3), '-o', 'Color', [0.9 0.3 0.6], 'LineWidth', 2 , 'MarkerFaceColor', [0.9 0.3 0.6], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,0+i) bstats(6,0+i)], ':', 'LineWidth', 3, 'Color', [0.9 0.3 0.6]);
end
%ci1 = plot(xaxes, bstats(10:12,23), 'o',  'Color', [0.9 0.3 0.6], 'MarkerFaceColor', [0.9 0.3 0.6]);
%ci2 = plot(xaxes, bstats(10:12,24), 'o',  'Color', [0.9 0.3 0.6], 'MarkerFaceColor', [0.9 0.3 0.6]);

xlim([0 1000]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
yticks([1.5:0.5:3.5]);
ylim([1 3.5]);
ax = gca;
ax.FontSize = 13; 
% yticks([1.5:0.5:3.5]);
% ylim([1 3.5]);
ylabel('Response onset (s)', 'fontsize', 15);

subplot(2,3,5)

xaxes = [100 300 900];
means = plot(xaxes, nbstats(1,4:6), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,3+i) nbstats(5,3+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(13:15,23), 'k.');
%ci2 = plot(xaxes, nbstats(13:15,24), 'k.');

means = plot(xaxes, bstats(1,4:6), '-o', 'Color', [0.4 0.9 0.9], 'LineWidth', 2, 'MarkerFaceColor', [0.4 0.9 0.9], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,3+i) bstats(6,3+i)], ':', 'LineWidth', 3, 'Color', [0.4 0.9 0.9]);
end
%ci1 = plot(xaxes, bstats(13:15,23), 'o', 'Color', [0.4 0.9 0.9], 'MarkerFaceColor', [0.4 0.9 0.9]);
%ci2 = plot(xaxes, bstats(13:15,24), 'o', 'Color', [0.4 0.9 0.9], 'MarkerFaceColor', [0.4 0.9 0.9]);

xlim([0 1000]);
xlabel('Stimulus duration (ms)', 'fontsize', 15);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
yticks([1.5:0.5:3.5]);
ylim([1 3.5]);
ax = gca;
ax.FontSize = 13; 


subplot(2,3,6)

means = plot(xaxes, nbstats(1,7:9), 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [nbstats(4,6+i) nbstats(5,6+i)],'LineWidth', 5, 'Color', [0.7 0.7 0.7]);
end
%ci1 = plot(xaxes, nbstats(16:18,23), 'k.');
%ci2 = plot(xaxes, nbstats(16:18,24), 'k.');

means = plot(xaxes, bstats(1,7:9), '-o', 'Color', [0.5 0.4 1], 'LineWidth', 2,'MarkerFaceColor', [0.5 0.4 1], 'MarkerSize', 7);
hold on
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,6+i) bstats(6,6+i)], ':', 'LineWidth', 3, 'Color', [0.5 0.4 1]);
end
%ci1 = plot(xaxes, bstats(16:18,23), 'o',  'Color', [0.5 0.4 1],'MarkerFaceColor', [0.5 0.4 1]);
%ci2 = plot(xaxes, bstats(16:18,24), 'o',  'Color', [0.5 0.4 1],'MarkerFaceColor', [0.5 0.4 1]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
xlim([0 1000]);

yticks([1.5:0.5:3.5]);
ylim([1 3.5]);

ax = gca;
ax.FontSize = 13; 

%% DURATION EFFECTS - PREDICTED VS OBSERVED
figure;

% Create standard impulse function (Friston)
imp_func = spm_hrf(0.001, [6 16 1 1 inf 0 32]);

subplot(2,3,1)
% Boxcars are defined in milliseconds. In this example, the first
% represents a stimulus of duration 1000 ms, amplitude 1. The second
% represents a stimulus of duration 100 ms, amplitude 10.
boxcar1 = [repmat(10,1,100) repmat(0,1,14900)];
boxcar2 = [repmat(10,1,300) repmat(0,1,14700)];
boxcar3 = [repmat(10,1,900) repmat(0,1,14100)];
hemo1 = conv(imp_func, boxcar1); %Perform the convolution of boxcar and standard impulse function
hemo2 = conv(imp_func, boxcar2);
hemo3 = conv(imp_func, boxcar3);
range = 1:15000; %Specify the timepoints of interest in hemo1 and hemo2 (in ms)
hemo1 = hemo1(range); %Trim the response to 24 seconds
hemo2 = hemo2(range);
hemo3 = hemo3(range);

%Plot the two hemodynamic responses
cols = {[0.9 0.3 0.6],[0.4 0.9 0.9],[0.5 0.4 1]};
plot(range/1000, hemo1, 'LineWidth', 2, 'Color', cols{1});
hold on;
plot(range/1000, hemo2, 'LineWidth', 2, 'Color', cols{2});
plot(range/1000, hemo3, 'LineWidth', 2, 'Color', cols{3});
yticks([0:0.4:2]);
ylim([-0.2 2]);
hold off
%Make everything look nicer
xlabel('Time (seconds) ');
ylabel('% signal change');
set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);
%Add a legend with information about the peak amplitude and response time
legend(['100ms'], ['300ms'], ['900ms']);
%Get the amplitude and time of the peak from each hemodynamic response
[peak1_amp, peak1_time] = max(hemo1);
[peak2_amp, peak2_time] = max(hemo2);
[peak3_amp, peak3_time] = max(hemo3);


subplot(2,3,2)

xaxes = [100 300 900];
amps = [peak1_amp, peak2_amp, peak3_amp];
plot(xaxes, amps, '-o', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerFaceColor', [0 0 0]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
xlim([0 1000]);
yticks([0:0.5:2]);
ylim([0 2]);
ylabel('Peak magnitude (%)');
xlabel('Stimulus duration (ms)');

set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);


subplot(2,3,3)
xaxes = [100 300 900];
times = [peak1_time, peak2_time, peak3_time] / 1000;
plot(xaxes, times, '-o', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerFaceColor', [0 0 0]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
xlim([0 1000]);
yticks([4:0.4:6.4]);
ylim([4 6.4]);
ylabel('Peak latency (s)');
xlabel('Stimulus duration (ms)');

set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);

subplot(2,3,4)

for l = 1:numel(lengths)
    for i = 3
        cols = {[0.9 0.3 0.6],[0.4 0.9 0.9],[0.5 0.4 1]};
        plot(0:0.625:15, eravecs{l,i}, 'Color', cols{l}, 'HandleVisibility','off')
        hold on
        TR = 0.625/1000;
        eventlength = 15*1000;
        t = 0:TR:(eventlength-1)*TR;
        t2 = 0:0.001:15-0.001;
        plot(t2, smoothhrfvec2{l,i}, 'Color', cols{l},'LineWidth', 2)
        %plot(t,hrfvecs{l,i},'Color', cols{l}, 'LineWidth', 2)
        ylim([-0.2 1.6]);
        yticks(0:0.4:1.6);
        xticks(0:2:16);
        xticklabels(xticklabs);
       
        legend('100ms', '300ms', '900ms', 'fontsize', 10);
        xlabel('Time (seconds) ');
        ylabel('% signal change');
        set(gca, 'fontunits', 'points', 'fontsize', 13);
        set(gca, 'fontunits', 'points', 'fontsize', 13);
        
    end
end

subplot(2,3,5)

bstats = csvread("C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out\auditory\MaxIntensity_BootStp_Amplitude_durationManipulation.csv",0,0);
xaxes = [100 300 900];
for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,0+i) bstats(6,0+i)], 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
hold on
end
hold on
means = plot(xaxes, bstats(1,1:3), '-o','Color', [0 0 0], 'LineWidth', 2, 'MarkerFaceColor', [0 0 0]);

xticks([100 300 900]);
xticklabels({'100', '300', '900'})
xlim([0 1000]);
yticks([0:0.5:2]);
ylim([0 2]);
ylabel('Peak magnitude (%)');
xlabel('Stimulus duration (ms)');
set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);


subplot(2,3,6)

bstats = csvread("C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out\auditory\MaxIntensity_BootStp_TimeToPeak_durationManipulation.csv",0,0);
xaxes = [100 300 900];

for i = 1:numel(xaxes)
bars = plot([xaxes(i) xaxes(i)], [bstats(5,0+i) bstats(6,0+i)], 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
hold on
end
hold on
means = plot(xaxes, bstats(1,1:3), '-o','Color', [0 0 0], 'LineWidth', 2, 'MarkerFaceColor', [0 0 0]);

%ci1 = plot(xaxes, bstats(16:18,7), 'o', 'Color', [0 0 0],'MarkerFaceColor', [0 0 0]);
%ci2 = plot(xaxes, bstats(16:18,8), 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
xticks([100 300 900]);
xticklabels({'100', '300', '900'})
xlim([0 1000]);
yticks([4:0.2:5]);
ylim([4 5]);
ylabel('Peak latency (s)');
xlabel('Stimulus duration (ms)');
set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);

hold off

%% INTENSITY EFFECTS - PREDICTED VS OBSERVED

figure;

subplot(2,3,1)
% Boxcars are defined in milliseconds. In this example, the first
% represents a stimulus of duration 1000 ms, amplitude 1. The second
% represents a stimulus of duration 100 ms, amplitude 10.
boxcar4 = [repmat(1,1,900) repmat(0,1,14100)];
boxcar5 = [repmat(3,1,900) repmat(0,1,14100)];
boxcar6 = [repmat(10,1,900) repmat(0,1,14100)];

hemo4 = conv(imp_func, boxcar4); %Perform the convolution of boxcar and standard impulse function
hemo5 = conv(imp_func, boxcar5);
hemo6 = conv(imp_func, boxcar6);
range = 1:15000; %Specify the timepoints of interest in hemo1 and hemo2 (in ms)
hemo4 = hemo4(range); %Trim the response to 24 seconds
hemo5 = hemo5(range);
hemo6 = hemo6(range);

%Plot the two hemodynamic responses
plot(range/1000, hemo4, 'LineWidth', 2, 'Color', [1 0 0]);
hold on;
plot(range/1000, hemo5, 'LineWidth', 2, 'Color', [0 1 0]);
plot(range/1000, hemo6, 'LineWidth', 2, 'Color', [0 0 1]);
yticks([0:0.4:2]);
ylim([-0.2 2]);


%Calculate plot area parameters to plot the boxcars nicely
a = axis;
maxAmp = 10;%max(max(boxcar1)-min(boxcar1), max(boxcar2)-min(boxcar2));

% %And now plot those boxcars
% plot((1:length(boxcar1))/1000, (boxcar1/maxAmp)*(a(4)/2.5)-a(4)/2, 'LineWidth', 1, 'Color', [1 0 0]);
% plot((1:length(boxcar2))/1000, (boxcar2/maxAmp)*(a(4)/2.5)-a(4)/3, 'LineWidth', 1, 'Color', [0 1 0]);
% plot((1:length(boxcar3))/1000, (boxcar3/maxAmp)*(a(4)/2.5)-a(4)/6, 'LineWidth', 1, 'Color', [0 0 0]);

%Make everything look nicer
set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);
xlabel('Time (seconds) ');
ylabel('% signal change');

%Add a legend with information about the peak amplitude and response time
legend(['10%'], ['30%'], ['90%']);


[peak4_amp, peak4_time] = max(hemo4);
[peak5_amp, peak5_time] = max(hemo5);
[peak6_amp, peak6_time] = max(hemo6);

subplot(2,3,2)

xaxes = [[log(0.01)/log(3)+1, log(0.1)/log(3)+1, log(1)/log(3)+1]];

amps = [peak4_amp, peak5_amp, peak6_amp];
plot(xaxes, amps, '-o', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerFaceColor', [0 0 0]);

xlim([-3.5 1.3]);
xlabv = [{'10', '30', '90'}];
xticks(xaxes);
xticklabels(xlabv);

% yticks([0:0.2:1.5]);
% ylim([0 1.5]);
yticks([0:0.5:2]);
ylim([0 2]);
ylabel('Peak magnitude (%)');
xlabel('Stimulus intensity (%)');

set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);


subplot(2,3,3)
xaxes = [[log(0.01)/log(3)+1, log(0.1)/log(3)+1, log(1)/log(3)+1]];

times = [peak4_time, peak5_time, peak6_time] / 1000;
plot(xaxes, times, '-o', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerFaceColor', [0 0 0]);

xlim([-3.5 1.3]);
xlabv = [{'10', '30', '90'}];
xticks(xaxes);
xticklabels(xlabv);

yticks([4:0.4:6.4]);
ylim([4 6.4]);
ylabel('Peak latency (s)');
xlabel('Stimulus intensity (%)');

set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);

subplot(2,3,4)

for l = 3
    for i = 1:numel(intensities)
        cols = ['r' 'g' 'b'];
        plot(0:0.625:15, eravecs{l,i}, cols(i),'HandleVisibility','off')
        hold on
        TR = 0.625/1000;
        eventlength = 15*1000;
        t = 0:TR:(eventlength-1)*TR;
        t2 = 0:0.001:15-0.001;
        plot(t2, smoothhrfvec2{l,i}, cols(i),'LineWidth', 2)
        %plot(t,hrfvecs{l,i}, cols(i), 'LineWidth', 2)
        ylim([-0.2 1.6]);
        yticks(0:0.4:1.6);
        xticks(0:2:16);
        xticklabels(xticklabs);
        legend('10%', '30%', '90%','fontsize', 10);
        xlabel('Time (seconds) ');
        ylabel('% signal change');
        set(gca, 'fontunits', 'points', 'fontsize', 13);
        set(gca, 'fontunits', 'points', 'fontsize', 13);
        
    end
end

subplot(2,3,5)
bstats = csvread("C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out\auditory\MaxDuration_BootStp_Amplitude_intensityManipulation.csv",0,0);

xaxes = [[log(0.01)/log(3)+1, log(0.1)/log(3)+1, log(1)/log(3)+1]];

for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [bstats(5,0+i) bstats(6,0+i)], 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
    hold on
end
hold on
means = plot(xaxes, bstats(1,1:3), '-o','Color', [0 0 0], 'LineWidth', 2, 'MarkerFaceColor', [0 0 0]);

xlim([-3.5 1.3]);
xlabv = [{'10', '30', '90'}];
xticks(xaxes);
xticklabels(xlabv);

% yticks([0:0.2:1.5]);
% ylim([0 1.5]);
yticks([0:0.5:2]);
ylim([0 2]);
ylabel('Peak magnitude (%)');
xlabel('Stimulus intensity (%)');
set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);


subplot(2,3,6)
bstats = csvread("C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out\auditory\MaxDuration_BootStp_TimeToPeak_intensityManipulation.csv",0,0);

xaxes = [[log(0.01)/log(3)+1, log(0.1)/log(3)+1, log(1)/log(3)+1]];

for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [bstats(5,0+i) bstats(6,0+i)], 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
    hold on
end
hold on
means = plot(xaxes, bstats(1,1:3), '-o','Color', [0 0 0], 'LineWidth', 2, 'MarkerFaceColor', [0 0 0]);

xlim([-3.5 1.3]);
xlabv = [{'10', '30', '90'}];
xticks(xaxes);
xticklabels(xlabv);

yticks([4:0.5:6]);
ylim([4 6]);

ylabel('Peak latency (s)');
xlabel('Stimulus intensity (%)');
set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);

hold off