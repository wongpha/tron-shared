%% TIMECOURSE EXTRACTION SCRIPT V2 [Alvin Wong, National University of Singapore]
%   This script extracts raw MR timecourse values averaged across voxels
%   within a specified ROI. In order to run this script, you must have FSL
%   installed on your system and have MATLAB launched from the Terminal. A
%   UNIX-based (e.g. MacOS) or Linux system is assumed.

% specify list of paths to input files (txt file)
%datadir = importdata('/Volumes/Passport/TRoN/TRoN_data/subs.txt');

subjectlist = importdata('/home/wongalvin/storage/SHARED_PROJECTS/BOLD6/PREPROC/subject_list.txt');

for subidx = 1:numel(subjectlist)
    currentsub = subjectlist(subidx, 1:7);
    datadir = dir(sprintf('%s', subjectlist(subidx, 8:end), '/vol/'), '*.nii.gz');
    
    for datidx = 1:numel(datadir) % step through each input
    %% ORGANIZATION OF SPECIFIED INPUTS    
        [currentdir,filename,~] = fileparts(datadir{datidx});
        [rootdir, currentdirname] = fileparts(currentdir);
        roidir = sprintf('%s', [rootdir, '/rois/']);
        %roidir = '/Volumes/Passport/TRoN/ROIs'; % paste the path to the folder where the ROI masks are saved here
        roilist = strsplit(strtrim(ls(roidir)))';
        inputfile = datadir{datidx};
        highres = [currentdir,'/norm_MNI152_1mm.nii.gz'];
        subjrunid = filename(1:14);
        intermediatename = filename(1:(end-4));

    %% UNCOMMENT BELOW LINES IF YOUR REGISTERED FUNCTIONAL DATA IS MISTAKENLY FLAGGED AS "SCANNER ANATOMICAL"
    %     fslorientargs = sprintf('%s %s', '-setsformcode 4', inputfile);
    %     fslorientcmd = sprintf('fslorient %s', fslorientargs);
    %     system(fslorientcmd);

    %% TIMECOURSE EXTRACTION LOOP [SEE PROGRESS IN COMMAND WINDOW]
        for roiidx = 1:numel(roilist)
            currentroi = sprintf('%s', [roidir, "/", roilist{roiidx}]);
            [roipath,roilabel,roiext] = fileparts(currentroi);
            roilabel=roilabel(1:(end-4));

            TCExtractFeedback = ['Currently extracting timecourses for', subjrunid, ' at ', roilabel];
            disp(TCExtractFeedback);

            outputfilename = sprintf('%s', [subjrunid, "_", roilabel]);
            TEargs = sprintf('%s %s %s %s %s', currentdir, inputfile, currentroi, outputfilename, subjrunid);
            TEcmd = sprintf('sh TCXfunc.sh %s', TEargs); % IMPT: TCXfunc.sh must be in same folder as this script
            system(TEcmd);
        end
    end

end