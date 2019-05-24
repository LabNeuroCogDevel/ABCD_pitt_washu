function [Framesretained, SubjectFD] = CalculateFD_267(subjectlist)
% Calculate FD across subjects
% Input: subjectlist = string containing path to txt file with sub list 

subjects = textread(subjectlist,'%s');

for s = 1:size(subjects,1)
    if s ~= 4
        subject = subjects{s};
        % Check # of runs 
        a = dir(['/data/nil-bluearc/GMT/ABCD/initial_ABCD_data/' subject '/MNINonLinear/Results/*rfMRI_REST*']);
        b = dir(['/data/nil-bluearc/GMT/ABCD/initial_ABCD_data/' subject '/MNINonLinear/Results/*rfMRI_REST_*']);
        n = numel(a) - numel(b);
        FD = [];
        for run = 1:n
            ThisFD = load(['/data/nil-bluearc/GMT/ABCD/initial_ABCD_data/' subject '/MNINonLinear/Results/rfMRI_REST' num2str(run) '/FNL_preproc_v2/FD.mat']);
            data = zeros(1+length(ThisFD.FD),1);
            data(1) = 0; data(2:end) = ThisFD.FD;
            FD = [FD; data];
        end
        SubjectFD(s,1) = mean(FD(2:end));
        FD(FD>=.2)=0; FD(FD<.2 & FD > 0) = 1; ThisFD = sum(FD);
        Framesretained(s,1) = ThisFD;
    else
        % Extract average time series for each ROI
        SubjectFD(s,1) = 1;
        Framesretained(s,1) = 1;  
    end

end

