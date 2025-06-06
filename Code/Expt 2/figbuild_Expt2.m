clear all
cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt2(2023)\GroupFit')
vstats = readmatrix('vstats.csv');
astats = readmatrix('astats.csv');
vgrouphrfdata = load('V_FittingMatrix.mat');
vgroupihrfdata = load('V_FittingMatrixI.mat');
agrouphrfdata = load('A_FittingMatrix.mat');
agroupihrfdata = load('A_FittingMatrixI.mat');
veravecs = load('V_ERAresults.mat');
aeravecs = load('A_ERAresults.mat');

%%
% This script takes the FittingMatrix cell from the TCParaFit group
% analysis outputs, which is saved as a MAT. It also takes the stats output
% from the PSTCParaFit script, which is saved as a CSV.

figure;
subplot(1,3,1);
    xaxes = [1.0 1.2 1.8];
    
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [vstats(4,3+i) vstats(5,3+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    
    means = plot(xaxes, vstats(1,4:6), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);

    xlim([0.8 2]);
    xticks([1.0 1.2 1.8]);
    xticklabels({'1.0','1.2','1.8'})
    yticks([0:0.2:2]);
    ylim([0 2]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak magnitude (% signal change)', 'fontsize', 14);
    %xlabel('Stimulus duration (s)', 'fontsize', 15);

    %title(sprintf('%s','Peak Magnitude'),'fontsize',16);

subplot(1,3,2);
    xaxes = [1.0 1.2 1.8];
   
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [vstats(4,0+i) vstats(5,0+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    
    means = plot(xaxes, vstats(1,1:3), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    
    xlim([0.8 2]);
    xticks([1.0 1.2 1.8]);
    xticklabels({'1.0','1.2','1.8'})
    yticks([5.2:0.1:6.4]);
    ylim([5.2 6.4]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak latency (s)', 'fontsize', 14);
    xlabel('Stimulus duration (s)', 'fontsize', 14);
    title(sprintf('%s','Visual Responses'),'fontsize',14);
    %title(sprintf('%s','Peak Latency'),'fontsize',16);

subplot(1,3,3);
    xaxes = [1.0 1.2 1.8];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [vstats(4,6+i) vstats(5,6+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, vstats(1,7:9), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([0.8 2]);
    xticks([1.0 1.2 1.8]);
    xticklabels({'1.0','1.2','1.8'})
    yticks([1.6:0.1:2.8]);
    ylim([1.6 2.8]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak onset (s)', 'fontsize', 14);
    %xlabel('Stimulus duration (s)', 'fontsize', 15);
    %title(sprintf('%s','Peak Onset'),'fontsize',16);



figure;
subplot(1,3,1);
    xaxes = [1.0 1.2 1.8];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [astats(4,3+i) astats(5,3+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, astats(1,4:6), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([0.8 2]);
    xticks([1.0 1.2 1.8]);
    xticklabels({'1.0','1.2','1.8'})
    yticks([0:0.2:2]);
    ylim([0 2]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak magnitude (% signal change)', 'fontsize', 14);
    %xlabel('Stimulus duration (s)', 'fontsize', 15);

    %title(sprintf('%s','Peak Magnitude'),'fontsize',16);

subplot(1,3,2);
    xaxes = [1.0 1.2 1.8];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [astats(4,0+i) astats(5,0+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, astats(1,1:3), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([0.8 2]);
    xticks([1.0 1.2 1.8]);
    xticklabels({'1.0','1.2','1.8'})
    yticks([4.0:0.1:5.2]);
    ylim([4.0 5.2]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak latency (s)', 'fontsize', 14);
    xlabel('Stimulus duration (s)', 'fontsize', 14);
    title(sprintf('%s','Auditory Responses'),'fontsize',14);
    %title(sprintf('%s','Peak Latency'),'fontsize',16);

subplot(1,3,3);
    xaxes = [1.0 1.2 1.8];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [astats(4,6+i) astats(5,6+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, astats(1,7:9), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([0.8 2]);
    xticks([1.0 1.2 1.8]);
    xticklabels({'1.0','1.2','1.8'})
    yticks([1.6:0.1:2.8]);
    ylim([1.6 2.8]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak onset (s)', 'fontsize', 14);
    %xlabel('Stimulus duration (s)', 'fontsize', 15);
    %title(sprintf('%s','Peak Onset'),'fontsize',16);

figure;
subplot(1,3,1);
    xaxes = [[log(0.1)/log(3)+1, log(1)/log(3)+1]];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [vstats(4,11+i) vstats(5,11+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, vstats(1,12:13), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([-1.4 1.3]);
    xlabv = [{'10', '100'}];
    xticks(xaxes);
    xticklabels(xlabv);

    yticks([0:0.2:2]);
    ylim([0 2]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak magnitude (% signal change)', 'fontsize', 14);
    
subplot(1,3,2);
    xaxes = [[log(0.1)/log(3)+1, log(1)/log(3)+1]];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [vstats(4,9+i) vstats(5,9+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, vstats(1,10:11), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([-1.4 1.3]);
    xlabv = [{'10', '100'}];
    xticks(xaxes);
    xticklabels(xlabv);
    
    yticks([5.2:0.1:6.4]);
    ylim([5.2 6.4]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak latency (s)', 'fontsize', 14);
    xlabel('Stimulus contrast (%)', 'fontsize', 14);
    title(sprintf('%s','Visual Responses'),'fontsize',14);
    
subplot(1,3,3);
    xaxes = [[log(0.1)/log(3)+1, log(1)/log(3)+1]];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [vstats(4,13+i) vstats(5,13+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, vstats(1,14:15), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([-1.4 1.3]);
    xlabv = [{'10', '100'}];
    xticks(xaxes);
    xticklabels(xlabv);
    yticks([1.6:0.1:2.8]);
    ylim([1.6 2.8]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak onset (s)', 'fontsize', 14);
    %xlabel('Stimulus duration (s)', 'fontsize', 15);
    %title(sprintf('%s','Peak Onset'),'fontsize',16);
    
    figure;
subplot(1,3,1);
    xaxes = [[log(0.1)/log(3)+1, log(1)/log(3)+1]];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [astats(4,11+i) astats(5,11+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, astats(1,12:13), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([-1.4 1.3]);
    xlabv = [{'30', '90'}];
    xticks(xaxes);
    xticklabels(xlabv);

    yticks([0:0.2:2]);
    ylim([0 2]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak magnitude (% signal change)', 'fontsize', 14);
    
subplot(1,3,2);
    xaxes = [[log(0.1)/log(3)+1, log(1)/log(3)+1]];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [astats(4,9+i) astats(5,9+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, astats(1,10:11), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([-1.4 1.3]);
    xlabv = [{'30', '90'}];
    xticks(xaxes);
    xticklabels(xlabv);
    
    yticks([4.0:0.1:5.2]);
    ylim([4.0 5.2]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak latency (s)', 'fontsize', 14);
    xlabel('Stimulus volume (%)', 'fontsize', 14);
    title(sprintf('%s','Auditory Responses'),'fontsize',14);
    
subplot(1,3,3);
    xaxes = [[log(0.1)/log(3)+1, log(1)/log(3)+1]];
    
    for i = 1:numel(xaxes)
    bars = plot([xaxes(i) xaxes(i)], [astats(4,13+i) astats(5,13+i)],'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
    hold on
    end
    hold on
    means = plot(xaxes, astats(1,14:15), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    hold on
    
    xlim([-1.4 1.3]);
    xlabv = [{'30', '90'}];
    xticks(xaxes);
    xticklabels(xlabv);
    yticks([1.6:0.1:2.8]);
    ylim([1.6 2.8]);
    % yticks([0:0.1:0.8]);
    % ylim([0 0.8]);
    ax = gca;
    ax.FontSize = 13; 
    ylabel('Peak onset (s)', 'fontsize', 14);
    %xlabel('Stimulus duration (s)', 'fontsize', 15);
    %title(sprintf('%s','Peak Onset'),'fontsize',16);
    %%
    
lengths = [1.0 1.2 1.8];
for l = 1:3
    for i = 2%:numel(intensities)
        
        hrfvecs{l,i}= vgrouphrfdata.FittingMatrix{l,i}.hrf;
        segment{l,i}= repmat(vgrouphrfdata.FittingMatrix{l,i}.hrf(end),1,1000);
        hrfvecs{l,i}= [hrfvecs{l,i} segment{l,i}];
%         erasegment{l,i} = veravecs.EventAvg{l,i}(end);
%         eravecs{l,i}= [veravecs.EventAvg{l,i}' erasegment{l,i}];
        ahrfvecs{l,i}= agrouphrfdata.FittingMatrix{l,i}.hrf;
        asegment{l,i} = repmat(agrouphrfdata.FittingMatrix{l,i}.hrf(end),1,1000);
        ahrfvecs{l,i}= [ahrfvecs{l,i} asegment{l,i}];
%         aerasegment{l,i} = aeravecs.EventAvg{l,i}(end);
%         auderavecs{l,i}= [aeravecs.EventAvg{l,i}' aerasegment{l,i}];
        
%         hrfvec2{l,i} = groupfittedhrfdata.FullHRFmatrix{l,i};
%         eravecs{l,i} = grouperadata.ERAmatrix{l,i};
%         smoothhrfvec2{l,i} = hrfvec2{l,i};
%         samplingRateIncrease = 1000;
%         newSamplePoints = linspace(0,15,15000);
%         smoothhrfvec2{l,i} = spline(0:0.625:15, smoothhrfvec2{l,i}, newSamplePoints);
    end
end

for l = 2
    for i = 1:2%:numel(intensities)
        
        hrfvec2{l,i}= vgroupihrfdata.FittingMatrix{l,i}.hrf;
        segment2{l,i} = repmat(vgroupihrfdata.FittingMatrix{l,i}.hrf(end),1,1000);
        hrfvec2{l,i}= [hrfvec2{l,i} segment2{l,i}];
%         erasegment2{l,i} = veravecs.EventAvg{l,i}(end);
%         eravecs2{l,i}= [veravecs.EventAvg{l,i}' erasegment{l,i}];
        ahrfvec2{l,i}= agroupihrfdata.FittingMatrix{l,i}.hrf;
        asegment2{l,i} = repmat(agroupihrfdata.FittingMatrix{l,i}.hrf(end),1,1000);
        ahrfvec2{l,i}= [ahrfvec2{l,i} asegment2{l,i}];
%         aerasegment2{l,i} = aeravecs.EventAvg{l,i}(end);
%         auderavecs2{l,i}= [aeravecs.EventAvg{l,i}' aerasegment{l,i}];

%         hrfvec2{l,i} = groupfittedhrfdata.FullHRFmatrix{l,i};
%         eravecs{l,i} = grouperadata.ERAmatrix{l,i};
%         smoothhrfvec2{l,i} = hrfvec2{l,i};
%         samplingRateIncrease = 1000;
%         newSamplePoints = linspace(0,15,15000);
%         smoothhrfvec2{l,i} = spline(0:0.625:15, smoothhrfvec2{l,i}, newSamplePoints);
    end
end
%%
figure;
subplot(1,2,1);
% Create standard impulse function (Friston)
imp_func = spm_hrf(0.001, [6 16 1 1 inf 0 32]);

% Boxcars are defined in milliseconds. In this example, the first
% represents a stimulus of duration 1000 ms, amplitude 1. The second
% represents a stimulus of duration 100 ms, amplitude 10.
boxcar1 = [repmat(10,1,1000) repmat(0,1,14000)];
boxcar2 = [repmat(10,1,1200) repmat(0,1,13800)];
boxcar3 = [repmat(10,1,1800) repmat(0,1,13200)];
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
yticks([0:0.5:3.5]);
ylim([-0.2 3.5]);

hold off
%Make everything look nicer
xlabel('Time (seconds) ');
ylabel('% signal change');
set(gca, 'fontunits', 'points', 'fontsize', 13);
set(gca, 'fontunits', 'points', 'fontsize', 13);
%Add a legend with information about the peak amplitude and response time
legend(['1000ms'], ['1200ms'], ['1800ms']);
title(sprintf('%s', ['Predicted']),'fontsize',16);
%Get the amplitude and time of the peak from each hemodynamic response
[peak1_amp, peak1_time] = max(hemo1);
[peak2_amp, peak2_time] = max(hemo2);
[peak3_amp, peak3_time] = max(hemo3);

subplot(1,2,2);
cols = {[0.5 0.5 0.5],[0 0 0]};
boxcar4 = [repmat(3,1,1200) repmat(0,1,13800)];
boxcar5 = [repmat(10,1,1200) repmat(0,1,13800)];

hemo4 = conv(imp_func, boxcar4); %Perform the convolution of boxcar and standard impulse function
hemo5 = conv(imp_func, boxcar5);

range = 1:15000; %Specify the timepoints of interest in hemo1 and hemo2 (in ms)
hemo4 = hemo4(range); %Trim the response to 24 seconds
hemo5 = hemo5(range);


%Plot the two hemodynamic responses
plot(range/1000, hemo4, 'LineWidth', 2, 'Color', cols{1});
hold on;
plot(range/1000, hemo5, 'LineWidth', 2, 'Color', cols{2});

yticks([0:0.5:3.5]);
ylim([-0.2 3.5]);


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
legend(['Low intensity'], ['High intensity']);
title(sprintf('%s', ['Predicted']),'fontsize',16);

[peak4_amp, peak4_time] = max(hemo4);
[peak5_amp, peak5_time] = max(hemo5);
%%
figure;

subplot(2,2,1);
    ybounds = [-0.4 1.6];
    for l = 1:3
        for i = 2
            cols = {[0.9 0.3 0.6],[0.4 0.9 0.9],[0.5 0.4 1]};
            %
            
            TR = 1/1000;
            eventlength = 15;
            t = -1:TR:eventlength-TR;
            t2 = 0:0.001:15-0.001;
            %plot(t2, smoothhrfvec2{l,i}, cols(i),'LineWidth', 2)
            plot(t, hrfvecs{l,i}, 'Color', cols{l}, 'LineWidth', 2)
            hold on
            plot(-1:1:15-1,veravecs.EventAvg{l,i}, 'Color', cols{l},'HandleVisibility','off')
            plot(14:15,repmat(veravecs.EventAvg{l,i}(end),2), 'Color', cols{l},'HandleVisibility','off')
            
            ylim(ybounds);
            yticks(ybounds(1):0.2:ybounds(2));
            xlim([-1 15]);
            xticks(-1:1:15);
            ax = gca;
            ax.FontSize = 12; 
            ylabel('% signal change', 'fontsize', 15);
            title(sprintf('%s', ['Visual']),'fontsize',16);
            legend('1s', '1.2s', '1.8s', 'fontsize', 15);
                        
        end
    end
    
subplot(2,2,2);
    ybounds = [-0.4 1.6];
    for l = 2
        for i = 1:2
            cols = {[0.5 0.5 0.5],[0 0 0]};
            
            TR = 1/1000;
            eventlength = 15;
            t = -1:TR:eventlength-TR;
            t2 = 0:0.001:15-0.001;
            %plot(t2, smoothhrfvec2{l,i}, cols(i),'LineWidth', 2)
            plot(t, hrfvec2{l,i}, 'Color', cols{i}, 'LineWidth', 2)
            hold on
            plot(-1:1:15-1,veravecs.EventAvg{l,i}, 'Color', cols{i},'HandleVisibility','off')
            plot(14:15,repmat(veravecs.EventAvg{l,i}(end),2), 'Color', cols{i},'HandleVisibility','off')
            
            ylim(ybounds);
            yticks(ybounds(1):0.2:ybounds(2));
            xlim([-1 15]);
            xticks(-1:1:15);
            ax = gca;
            ax.FontSize = 12; 
            %ylabel('% signal change', 'fontsize', 15);
            
            title(sprintf('%s', ['Visual']),'fontsize',16);
            legend('10%', '100%', 'fontsize', 15);
        end
    end
    
subplot(2,2,3);
    ybounds = [-0.4 1.6];
    for l = 1:3
        for i = 2
            cols = {[0.9 0.3 0.6],[0.4 0.9 0.9],[0.5 0.4 1]};
            
            TR = 1/1000;
            eventlength = 15;
            t = -1:TR:eventlength-TR;
            t2 = 0:0.001:15-0.001;
            %plot(t2, smoothhrfvec2{l,i}, cols(i),'LineWidth', 2)
            plot(t, ahrfvecs{l,i}, 'Color', cols{l}, 'LineWidth', 2)
            hold on
            plot(-1:1:15-1,aeravecs.EventAvg{l,i}, 'Color', cols{l},'HandleVisibility','off')
            plot(14:15,repmat(aeravecs.EventAvg{l,i}(end),2), 'Color', cols{l},'HandleVisibility','off')
            
            ylim(ybounds);
            yticks(ybounds(1):0.2:ybounds(2));
            xlim([-1 15]);
            xticks(-1:1:15);
            ax = gca;
            ax.FontSize = 12; 
            ylabel('% signal change', 'fontsize', 15);
            xlabel('Time since event onset (s)', 'fontsize', 15);
            title(sprintf('%s', ['Auditory']),'fontsize',16);
            legend('1s', '1.2s', '1.8s', 'fontsize', 15);
        end
    end
    
subplot(2,2,4);
    ybounds = [-0.4 1.6];
    for l = 2
        for i = 1:2
            cols = {[0.5 0.5 0.5],[0 0 0]};
            TR = 1/1000;
            eventlength = 15;
            t = -1:TR:eventlength-TR;
            t2 = 0:0.001:15-0.001;
            %plot(t2, smoothhrfvec2{l,i}, cols(i),'LineWidth', 2)
            plot(t, ahrfvec2{l,i}, 'Color', cols{i}, 'LineWidth', 2)
            hold on
            
            plot(-1:1:15-1,aeravecs.EventAvg{l,i}, 'Color', cols{i},'HandleVisibility','off')
            plot(14:15,repmat(aeravecs.EventAvg{l,i}(end),2), 'Color', cols{i},'HandleVisibility','off')
            
            ylim(ybounds);
            yticks(ybounds(1):0.2:ybounds(2));
            xlim([-1 15]);
            xticks(-1:1:15);
            ax = gca;
            ax.FontSize = 12; 
            %ylabel('% signal change', 'fontsize', 15);
            xlabel('Time since event onset (s)', 'fontsize', 15);
            title(sprintf('%s', ['Auditory']),'fontsize',16);
            legend('30%', '90%', 'fontsize', 15);
        end
    end