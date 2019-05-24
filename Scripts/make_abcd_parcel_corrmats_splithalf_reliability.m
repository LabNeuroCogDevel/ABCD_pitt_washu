function make_abcd_parcel_corrmats_splithalf_reliability(subjects,numpools,inputdir)
% Loop through a subject list, import time series, extract parcel time series,
% frame censor, make correlation matrices 
% Parcels are from Gordon & Laumann et al 2017

% Inputs
%   subjects = cell - each contains a unique subject string 
%   inputdir = string - input directory 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts'))

Corrmat = zeros(333,333,2,size(subjects,1));
Reliability = nan(size(subjects,1),1);
%MissingSubIdx = zeros(size(subjects,1),1);
%Badtmaskidx = zeros(size(subjects,1),1);
%SubjectFD = zeros(size(subjects,1),1);
%SubjectFD_Incframes = zeros(size(subjects,1),1);
%SubjectFrameTotal = zeros(size(subjects,1),1);

%q = parpool(numpools);

for s = 1:size(subjects,1)
    fileexists = 0;
    subject = subjects{s};

    inputfile = [inputdir 'sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/task-rest_DCANBOLDProc_v4.0.0_Gordon.ptseries.nii'];
    dir_badtask = ['/data/Daenerys/ABCD/data/2k_bad_task/sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/task-rest_DCANBOLDProc_v4.0.0_Gordon.ptseries.nii'];
    if exist(inputfile)
        path = inputdir;
        fileexists = 1;
    elseif exist(dir_badtask)
        path = '/data/Daenerys/ABCD/data/2k_bad_task/';
        fileexists = 1;
    else
        disp(['Subject ' num2str(s) ' is missing!'])
        MissingSubIdx(s,1) = 1;
    end
    
    if fileexists
        disp(['Subject ' num2str(s)]) 
        ABCD = ft_read_cifti_mod([path 'sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/task-rest_DCANBOLDProc_v4.0.0_Gordon.ptseries.nii']);
        CortData = ABCD.data; 
        % Temp masks
        % Check # of runs 

        n = numel(dir([path 'sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/*task-rest0*']));
        FD = [];
        for run = 1:n
            if exist([path 'sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/task-rest0' num2str(run) '/DCANBOLDProc_v4.0.0/FD.mat'])
                ThisFD = load([path 'sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/task-rest0' num2str(run) '/DCANBOLDProc_v4.0.0/FD.mat']);
                data = zeros(1+length(ThisFD.FD),1);
                data(1) = 0; data(2:end) = ThisFD.FD;
                FD = [FD; data];
            end
        end
        %SubjectFD(s,1) = mean(FD(2:end));
        %tmp = FD(2:end);
        %SubjectFD_Incframes(s,1) = mean(tmp(tmp<.2)); clear tmp
        FD(FD>=.2)=0; FD(FD<.2 & FD > 0) = 1; ThisFD = sum(FD);
        %SubjectFrameTotal(s,1) = ThisFD;
        % apply temporal mask. 
        FD = repmat(FD,1,size(CortData,1))';
        if size(FD,2) == size(CortData,2)
           CortData=CortData.*FD;
           CortData(:,~any(CortData,1)) = [];    
        else
           disp(['Subject ' num2str(s) 'has tmask mismatch'])
           %Badtmaskidx(s,1) = 1;
        end

        FH = FisherTransform(corr(CortData(:,1:(floor(size(CortData,2)/2)))'));
        uidx = find(triu(FH,1));
        SH = FisherTransform(corr(CortData(:,(1+ceil(size(CortData,2)/2)):size(CortData,2))'));
        Reliability(s,1) = corr(FH(uidx),SH(uidx));
        Corrmat(:,:,1,s) = FH;
        Corrmat(:,:,2,s) = SH;
        
        disp(['Reliability for subject ' num2str(s) ' = ' num2str(Reliability(s,1))])
        if s > 1
            disp(['Mean Reliability = ' num2str(nanmean(Reliability(1:s,1)))])
        end
        clear FH SH CortData FD ABCD
    end

end

disp('           **************   Saving   ***********            ')
save('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SplitHalves_SubjectParcelCorrMats.mat','Corrmat','-v7.3')
save('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SplitHalfReliability_rest.mat','Reliability')

end

