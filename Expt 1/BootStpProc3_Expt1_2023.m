clearvars();
cd('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\Data_Pub');
eventtype = "A";
filedir = dir('A1*.mat');

for fileidx = 1:numel(filedir)
    clearvars('-except', 'eventtype', 'filedir', 'fileidx', 'SubjAvgTC', 'PerSubjBootStpInput');
    filename = sprintf('%s', [filedir(fileidx).folder, '/', filedir(fileidx).name]);
    load(filename);
    subjname = [filedir(fileidx).name(4:10)];
    
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
    elseif eventtype == 'A'
        conddefinition = conddef{1};
        eventlengths = [0.1 0.3 0.9];
        eventintensities = [0.1 0.3 1];
        roiname = 'A1';
    end            

    PerSubjBootStpInput{fileidx} = ProtoMatrix;
    
end

%% Uncomment to generate new bootstrap sample
%[bootstat, bootsam_raw] = bootstrp(10000, [], [1:numel(filedir)]);
%save("bootstrap_sample", 'bootsam_raw');
load('C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\Data_Pub\bootstrap_sample.mat');
bootsam = bootsam_raw;
%bootsam = reshape(bootsam_raw, 15, 100, []);


%for h = 1:size(bootsam,3) % pages
    
    %clearvars('OutputMatrix');
    InputMatrix{3, 3} = [];
    OutputMatrix{3, 3} = [];
    
    for x = 1:size(OutputMatrix,1) % init output
        for y = 1:size(OutputMatrix,2)
            OutputMatrix{x,y} = OutputStructInit();
        end
    end
    
    for i = 1:size(bootsam,2) % iterations
        for j = 1:size(bootsam,1) % members
            %bootstrap_sample{j} = PerSubjBootStpInput{bootsam(j,i,h)};
            
            bootstrap_sample{j} = PerSubjBootStpInput{bootsam(j,i)};
            
            for k=1:size(bootstrap_sample{j},1) % conditions
                for l=1:size(bootstrap_sample{j},2)
                    bootstrap_sample{j}{k,l} = reshape(bootstrap_sample{j}{k,l}, size(bootstrap_sample{j}{k,l},1), size(bootstrap_sample{j}{k,l},2)*size(bootstrap_sample{j}{k,l},3));
                    BaselineVector{j}{k,l} = bootstrap_sample{j}{k,l}(13, :);
                    for m=1:size(bootstrap_sample{j}{k,l}, 2)
                        bootstrap_sample{j}{k,l}(:,m) = ((bootstrap_sample{j}{k,l}(:,m) - BaselineVector{j}{k,l}(1,m)) / BaselineVector{j}{k,l}(1,m)) * 100;
                    end
                    bootstrap_sample{j}{k,l} = nanmean(bootstrap_sample{j}{k,l},2); %% comment out if you don't want timecourses to be averaged within subjects
                    InputMatrix{k,l}(:,j) = bootstrap_sample{j}{k,l};
                end
            end
        end

        for lengths = 1:size(InputMatrix,1)
            for ints = 1:size(InputMatrix,2)
                [FittingMatrix{lengths, ints}] = GammaFit(InputMatrix{lengths, ints}(10:34,:)', 0.625, 15); % event onsets at index 10
                FittingMatrix{lengths, ints}.StimDur = eventlengths(lengths);
                FittingMatrix{lengths, ints}.StimInt = eventintensities(ints);
                FittingMatrix{lengths, ints}.Iter = i;
                FittingMatrix{lengths, ints}.BootSam = bootsam(:,i);
                %FittingMatrix{lengths, ints}.Iter = [h,i];
                disp(['Iteration: ', int2str(FittingMatrix{lengths, ints}.Iter)]);
                OutputMatrix{lengths, ints}(size(OutputMatrix{lengths, ints},2) + 1) = FittingMatrix{lengths,ints};
            end
        end


    end

    for x = 1:size(OutputMatrix,1) % clean up output
        for y = 1:size(OutputMatrix,2)
            OutputMatrix{x,y} = OutputMatrix{x,y}(2:end);
        end
    end
    
    outputdir = "C:\Users\alvinw\iCloudDrive\Analysis\TronAnalysis\Expt1 [2021]\PSubj2\boot_out";
    save(sprintf('%s', outputdir, filesep, 'OutputMatrix'), "OutputMatrix");
    %save(sprintf('%s', outputdir, filesep, 'OutputMatrix_cluster', int2str(h)), "OutputMatrix");
    
%end