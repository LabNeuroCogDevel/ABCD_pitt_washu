function make_abcd_vertex_corrmats(subjects,inputdir,outname)
% Loop through a subject list, import time series, extract parcel time series,
% frame censor, make correlation matrices 
% Parcels are from Gordon & Laumann et al 2017

% Inputs
%   subjects = cell - each contains a unique subject string 
%   inputdir = string - input directory 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts'))


MissingSubIdx = zeros(size(subjects,1),1);
Badtmaskidx = zeros(size(subjects,1),1);
SubjectFD = zeros(size(subjects,1),1);
SubjectFrameTotal = zeros(size(subjects,1),1);

%pool = parpool(numpools);
for s = 1:size(subjects,1)
    fileexists = 0;
    subject = subjects{s};

    inputfile = [inputdir 'sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/task-rest_DCANBOLDProc_v4.0.0_Atlas.dtseries.nii'];
    dir_badtask = ['/data/Daenerys/ABCD/data/2k_bad_task/sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/task-rest_DCANBOLDProc_v4.0.0_Atlas.dtseries.nii'];
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
        
        ABCD = ft_read_cifti_mod([path 'sub-NDAR' subject '/ses-baselineYear1Arm1/files/MNINonLinear/Results/task-rest_DCANBOLDProc_v4.0.0_Atlas.dtseries.nii']);
        %brainstructind = ABCD.brainstructure(ABCD.brainstructure>0);
        CortData = ABCD.data(1:59412,:); 

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
        SubjectFD(s,1) = mean(FD(2:end));
        FD(FD>=.2)=0; FD(FD<.2 & FD > 0) = 1; ThisFD = sum(FD);
        SubjectFrameTotal(s,1) = ThisFD;
        % apply temporal mask. 
        FD = repmat(FD,1,size(CortData,1))';
        if size(FD,2) == size(CortData,2)
            CortData=CortData.*FD;
            CortData(:,~any(CortData,1)) = [];
        else
            disp(['Subject ' num2str(s) 'has tmask mismatch'])
            Badtmaskidx(s,1) = 1;
        end
        clear FD
        if s == 1
            Corrmat = zeros(size(CortData,1));
            Corrmat = FisherTransform(corr(CortData'));
            clear CortData
        else
            data = FisherTransform(corr(CortData'));
            Corrmat = (((s-1).*Corrmat) + data) ./ (s);
            clear data
        end
        
        
    end

    if floor(s/100)==s/100 % save every 100 subjects 
        disp('           **************   Saving   ***********          ')
        %save('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD.mat','SubjectFD')
        %save('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFrameTotal.mat','SubjectFrameTotal')
        %save('/data/nil-bluearc/GMT/Scott/ABCD/FD/Badtmaskidx.mat','Badtmaskidx')
        save('/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/MissingSubIdx.mat','MissingSubIdx')  
        save('/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/MeanCorrmat.mat','Corrmat','-v7.3')
    end

end

disp('           **************   Saving   ***********            ')
%save('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD.mat','SubjectFD')
%save('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFrameTotal.mat','SubjectFrameTotal')
%save('/data/nil-bluearc/GMT/Scott/ABCD/FD/Badtmaskidx.mat','Badtmaskidx')
save('/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/MissingSubIdx.mat','MissingSubIdx')
save('/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/MeanCorrmat.mat','Corrmat','-v7.3')

% Save dconn
output = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/ABCDtest.dconn.nii');
output.data(1:59412,1:59412) = single(Corrmat);
ft_write_cifti_mod(outname,output)


end