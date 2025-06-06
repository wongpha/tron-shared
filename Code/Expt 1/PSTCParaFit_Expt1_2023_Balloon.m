clearvars();

cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\Data'); % Path to data matrices output by PerSubjProc
eventtype = "V"; % Set ROI type (A or V)
filedir = dir('V1*.mat'); % Set filename search criteria
varNames = {'TTP1', 'TTP2', 'TTP3', 'AMP1', 'AMP2', 'AMP3', 'ONS1', 'ONS2', 'ONS3'};
varTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
stattable{1} = table('Size', [numel(filedir), 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{2} = table('Size', [numel(filedir), 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{3} = table('Size', [numel(filedir), 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);

stattable{4} = table('Size', [numel(filedir), 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{5} = table('Size', [numel(filedir), 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
stattable{6} = table('Size', [numel(filedir), 3*3], 'VariableTypes', varTypes, 'VariableNames', varNames);
global signal
global covs

for fileidx = 1:numel(filedir)
    clearvars('-except','signal', 'covs', 'varNames', 'varTypes', 'stattable',... 
        'eventtype', 'filedir', 'customvolsperevent', 'volsperevent', 'fileidx', 'AllStats', 'SubjFittedHRFs', 'SubjFittedFullHRFs', 'SubjERAs');
    filename = sprintf('%s', [filedir(fileidx).folder, '/', filedir(fileidx).name]);
    load(filename);
    subjname = [filedir(fileidx).name(4:10)];

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
        outputdir = "C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\out\vis_balloon"
        
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
        outputdir = "C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\out\aud_balloon"

        for datacnt = 1:numel(CurrentA)
        AOnset(:, datacnt) = CurrentA(datacnt).AuditoryEvs(1:end);
        ALength(:,datacnt) = CurrentA(datacnt).AuditoryLengths;
        AInt(:,datacnt) = CurrentA(datacnt).AuditoryInts;
        end
        
        for datacnt = 1:numel(CurrentA)
        for evcnt = 1:size(AOnset(:,datacnt),1)
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
    
    for conC = 1:numel(ProtoMatrix) % baseline computation
        BaselineVector{conC} = ProtoMatrix{conC}(13, :);
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
    

    % Per-condition balloon function fitting for each subject

    for lengths = 1:size(ProtoMatrix,1)
        for ints = 1:size(ProtoMatrix,2)
            
            input = ProtoMatrix{lengths,ints}(10:34,:); % 20s timecourses
            trials = size(ProtoMatrix{lengths,ints}, 2); % total number of events
            totalcourselength = numel(ProtoMatrix{lengths,ints}); % calculate total number of volumes going in
            timecourselength = 15; % length of each timecourse (s)
            numvols = 25; % number of volumes per event
            
            % simulate SPM structs
            SPM.xY.RT = 0.625; % repetition time (s)
            SPM.nscan = trials*numvols; % total no. of volumes (num of events * num of vols per event)
            SPM.xBF.T = 6.25; % no. of subdivisions of TR
            SPM.xBF.name = 'hrf';
            SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T; % length of each subdivision (s)
            SPM.xBF.UNITS = 'secs';
            SPM.Sess.U.name = {sprintf('%s', ['vis', '_','l', num2str(lengths), 'i', num2str(ints)])}; % name of condition
            SPM.Sess.U.u = []; % this is required
            SPM.Sess.U.ons = [3:timecourselength:totalcourselength*SPM.xY.RT]; % event onsets in UNITS
            SPM.Sess.U.dur = 0; %eventlengths(lengths); %0; % set 0 for event-related designs (counter-intuitive I know; see section 8.2.1 of SPM manual)
            SPM.Sess.U.orth = 1; % orthogonality of EV
            SPM.Sess.U = spa_get_ons(SPM,1); % builds condition specification from info thus far

            s = length(SPM.Sess);
            Sess = SPM.Sess(s);
            U.dt = Sess.U(1).dt;
            u    = length(Sess.U);
            U.name = SPM.Sess.U.name;
            U.u    = Sess.U(1).u(33:end,1);

            TE    = 0.04;

            Y.y = reshape(input, [trials*numvols, 1]); % converts protomatrix into one long timecourse
            Y.dt = 0.625; %SPM.xY.RT;

            % Model specification: m input; 4 states; 1 output; m + 6 parameters
            %--------------------------------------------------------------------------
            % u(m) - mth stimulus function     (u)
            %
            % x(1) - vascular signal           log(s)
            % x(2) - rCBF                      log(f)
            % x(3) - venous volume             log(v)
            % x(4) - deoyxHb                   log(q)
            %
            % y(1) - BOLD                      (y)
            %
            % P(1)       - signal decay               d(ds/dt)/ds)      half-life = log(2)/P(1) ~ 1sec
            % P(2)       - autoregulation             d(ds/dt)/df)      2*pi*sqrt(1/P(1)) ~ 10 sec
            % P(3)       - transit time               (t0)              ~ 1 sec
            % P(4)       - exponent for Fout(v)       (alpha)           c.f. Grubb's exponent (~ 0.38)
            % P(5)       - resting oxygen extraction  (E0)              ~ range 20 - 50%
            % P(6)       - ratio of intra- to extra-  (epsilon)         ~ range 0.5 - 2
            %              vascular components   
            %              of the gradient echo signal   

            % P(6 + 1:m) - input efficacies - d(ds/dt)/du)  ~0.3 per event
            %--------------------------------------------------------------------------

            % priors (3 modes of hemodynamic variation)
            %--------------------------------------------------------------------------
            m       = size(U.u,2);
            [pE,pC] = spa_hdm_priors(m,3);

            % model
            %--------------------------------------------------------------------------
            M.f     = 'spa_fx_hdm';
            M.g     = 'spa_gx_hdm';
            M.x     = [0 0 0 0]'; 
            M.pE    = pE;    
            M.pC    = pC;
            M.m     = m;
            M.n     = 4;
            M.l     = 1;
            M.N     = 1500; % upscaling of TR for smooth model fits
            M.dt    = timecourselength/M.N; % event lengths (s) / M.N
            M.TE    = TE;
            
            % nonlinear system identification
            %--------------------------------------------------------------------------
            [Ep,Cp,Ce,K0,K1,K2,M0,M1,L1,L2,F] = spa_nlsi(M,U,Y);
            
            t       = [1:M.N]*M.dt;
            Fhdm    = spm_figure;
            set(Fhdm,'name','Hemodynamic Modeling')

            % display input parameters
            %--------------------------------------------------------------------------
            subplot(2,2,1)
            P     = Ep(7:end);
            C     = diag(Cp(7:end,7:end));
            [i,j] = max(abs(P));
            spm_barh(P,C)
            axis square
            title({'stimulus efficacy'; 'with 90% confidence intervals'},'FontSize',10)
            set(gca,'Ytick',[1:m],'YTickLabel',U.name,'FontSize',8)
            str = {};
            for i = 1:m
                str{end + 1} = U.name{i};
                str{end + 1} = sprintf('mean = %0.2f',P(i));
                str{end + 1} = '';
            end
            set(gca,'Ytick',[1:m*3]/3 + 1/2,'YTickLabel',str)
            xlabel('relative efficacy per event/sec')

            % display hemodynamic parameters
            %--------------------------------------------------------------------------
            subplot(2,2,3)
            P     = Ep(1:6);
            pE    = pE(1:6);
            C     = diag(Cp(1:6,1:6));
            spm_barh(P,C,pE)
            title({ 'hemodynamic parameters'},'FontSize',10)
            set(gca,'Ytick',[1:18]/3 + 1/2)
            set(gca,'YTickLabel',{  'SIGNAL decay',...
                        sprintf('%0.2f per sec',P(1)),'',...
                        'FEEDBACK',...
                        sprintf('%0.2f per sec',P(2)),'',...
                        'TRANSIT TIME',...
                        sprintf('%0.2f seconds',P(3)),'',...
                        'EXPONENT',...
                        sprintf('%0.2f',P(4)),'',...
                        'EXTRACTION',...
                        sprintf('%0.0f %s',P(5)*100,'%'),'',...
                        'log SIGNAL RATIO',...
                        sprintf('%0.2f %s',P(6),'%'),''},'FontSize',8)


            % get display state kernels (i.e. state dynamics) 
            %==========================================================================

            % Volterra kernels of states
            %--------------------------------------------------------------------------
            [H0,H1] = spa_kernels(M0,M1,M.N,M.dt);

            subplot(3,2,2)
            plot(t,exp(H1(:,:,j)))
            axis square
            title({['1st order kernels for ' U.name{j}];...
                'state variables'},'FontSize',9)
            ylabel('normalized values')
            legend({'s','f','v','q'}, 'Location','Best');
            grid on


            % display output kernels (i.e. BOLD response) 
            %--------------------------------------------------------------------------
            subplot(3,2,4)
            plot(t,K1(:,:,j))
            axis square
            title({'1st order kernel';...
                'output: BOLD'},'FontSize',9)
            ylabel('normalized flow signal')
            grid on

            subplot(3,2,6)
            imagesc(t,t,K2(:,:,1,j,j))
            axis square
            title({'2nd order kernel';...
                'output: BOLD'},'FontSize',9)
            xlabel({'time \{seconds\} for'; U.name{j}})
            grid on
            %%
            H2 = exp(H1);
            Hsave{lengths,ints}(:,:,fileidx) = H2;
            Hmean{lengths,ints} = nanmean(Hsave{lengths,ints},3);
            FPcovariance{lengths,ints}(1:numel(covs),fileidx) = covs;
            p = 0.4;
            v0 = 0.02;
            nu0   = 40.3;
            r0 = 25;
            E0 = P(5);
            epsi = exp(P(6));
            k1 = 6.7; %4.3.*nu0.*E0.*TE;
            k2 = 2.73; %epsi.*r0.*E0.*TE; 
            k3    = 1-epsi;
            modelvoltkernel1{lengths,ints}(:,fileidx) = 100/p*v0*(k1.*(1-Hmean{lengths,ints}(:,4))+k2.*(1-Hmean{lengths,ints}(:,4)./Hmean{lengths,ints}(:,3)+k3.*(1-Hmean{lengths,ints}(:,3))));
            modelvoltkernel2{lengths,ints}(:,fileidx) = K1;
            fullEMPTC{lengths,ints}(1:numel(signal)) = signal;
            close
        end
    end
    for lengths = 1:size(ProtoMatrix,1)
        for ints = 1:size(ProtoMatrix,2)
            Volkernel1{lengths, ints} = nanmean(modelvoltkernel1{lengths, ints},2);
            Volkernel2{lengths, ints} = nanmean(modelvoltkernel2{lengths, ints},2);
            rsEMPTC{lengths,ints} = reshape(fullEMPTC{lengths,ints}, numvols, numel(fullEMPTC{lengths,ints})/numvols);
            EMPredict{lengths,ints} = nanmean(rsEMPTC{lengths, ints}(:,2:end),2);
            EMPredictSeg{lengths,ints} = EMPredict{lengths,ints};
            EMPredictSmooth{lengths,ints}.hrf = resample(EMPredict{lengths,ints}, 1000,1);
            [EMPredictSmooth{lengths,ints}.MaxValue, EMPredictSmooth{lengths,ints}.MaxIndex] = max(EMPredictSmooth{lengths,ints}.hrf);
            EMPredictSmooth{lengths,ints}.MaxIndexInSeconds = EMPredictSmooth{lengths,ints}.MaxIndex / 1000 * 0.625;
            
            EMPredictSmooth{lengths, ints}.HalfMax = (EMPredictSmooth{lengths,ints}.MaxValue-0)/2 + 0;
            EMPredictSmooth{lengths,ints}.OnsetMax = (EMPredictSmooth{lengths,ints}.MaxValue-0)/10 + 0;
            [EMPredictSmooth{lengths,ints}.MinValue, EMPredictSmooth{lengths,ints}.MinIndex] = min(EMPredictSmooth{lengths,ints}.hrf);

            for i = EMPredictSmooth{lengths,ints}.MaxIndex:-1:1
                EMPredictSmooth{lengths,ints}.HMPoint1 = i;
                if EMPredictSmooth{lengths,ints}.hrf(i) < EMPredictSmooth{lengths,ints}.HalfMax
                    break
                end
            end

            for i = EMPredictSmooth{lengths,ints}.MaxIndex:24000
                EMPredictSmooth{lengths,ints}.HMPoint2 = i;
                if EMPredictSmooth{lengths,ints}.hrf(i) < EMPredictSmooth{lengths,ints}.HalfMax
                    break
                end
            end

            for i = EMPredictSmooth{lengths,ints}.MaxIndex:-1:1
                EMPredictSmooth{lengths,ints}.RisePoint = i;
                if EMPredictSmooth{lengths,ints}.hrf(i) < EMPredictSmooth{lengths,ints}.OnsetMax
                    break
                end
            end

            EMPredictSmooth{lengths,ints}.FWHM = EMPredictSmooth{lengths,ints}.HMPoint2 - EMPredictSmooth{lengths,ints}.HMPoint1;
            EMPredictSmooth{lengths,ints}.OnsetTimeInSeconds = EMPredictSmooth{lengths,ints}.RisePoint / 1000 * 0.625;
            EMPredictSmooth{lengths,ints}.DispersionInSeconds = EMPredictSmooth{lengths,ints}.FWHM / 1000 * 0.625;
            
            PeakTimeMatrix{lengths,ints} = EMPredictSmooth{lengths,ints}.MaxIndexInSeconds;
            PeakAmplitudeMatrix{lengths,ints} = EMPredictSmooth{lengths,ints}.MaxValue;
            FWHMMatrix{lengths,ints} = EMPredictSmooth{lengths,ints}.DispersionInSeconds;
            PeakOnsetMatrix{lengths, ints} = EMPredictSmooth{lengths,ints}.OnsetTimeInSeconds;
            HRFMatrix{lengths,ints} = EMPredictSmooth{lengths,ints}.hrf;
        end
    end
    FittingMatrix = [EMPredictSmooth, HRFMatrix];

    cd(outputdir)
    save(sprintf('%s', [eventtype, '_', subjname, '_FittingMatrix_Balloon', '.mat']),'FittingMatrix');

    AllStats{fileidx,1} = [{PeakTimeMatrix}, {PeakAmplitudeMatrix}, {PeakERA}, {TimeToPeakERA}, {FWHMMatrix}, {PeakOnsetMatrix}];
    
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
