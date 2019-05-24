% PFactor analysis for 267 subjects 
% Path to behavioral data
Behavioraldata = '/data/nil-bluearc/GMT/Scott/ABCD/Pfactor_onefactor_allsubjects_4k.txt';
% Read in subject list for brain & behavior 
big_subject_list = '/data/nil-bluearc/GMT/Scott/ABCD/ABCD_subjectlist_pfactor_4k.txt';
small_subject_list = '/data/nil-bluearc/GMT/Scott/ABCD/ABCD_subjectlist_267.txt';
SubjectMatchIdx = MatchSubjectList(big_subject_list,small_subject_list);

% Load p-factor data (one factor only)
Data = textread(Behavioraldata,'%s');
for i = 1:size(Data,1)
    Pscores(i,1) = str2num(Data{i});
end

PscoresSublist = Pscores(SubjectMatchIdx==1);

% Read in 267 individual corr mats 
load('/data/nil-bluearc/GMT/Scott/ABCD_DCN_CognitionPaper/SubjectCorrMats.mat')

%% Calculate FD and kick out high motion (>.2) subs 
[Framesretained, SubjectFD] = CalculateFD_267(subjectlist);


%%
LowFDidx = find(SubjectFD<=.2 & Framesretained>(480/.8)); % mean FD <= 0.2 and 8 mins of data
MatchingCorrMats = Corrmat(:,:,LowFDidx);
MatchingPfactorscores = PscoresSublist(LowFDidx);
BetaMat = zeros(size(MatchingCorrMats,1),size(MatchingCorrMats,2));
TMat = zeros(size(MatchingCorrMats,1),size(MatchingCorrMats,2));
rMat = zeros(size(MatchingCorrMats,1),size(MatchingCorrMats,2));
for i = 1:(size(MatchingCorrMats,1)-1)
    disp(['On ROI ' num2str(i)])
    
    for j = i+1:size(MatchingCorrMats,1)

        a = LinearModel.fit(squeeze(MatchingCorrMats(i,j,:)),MatchingPfactorscores);
        BetaMat(i,j) = table2array(a.Coefficients(2,1)); % beta
        TMat(i,j) = table2array(a.Coefficients(2,3)); % t stat
        rMat(i,j) = corr(squeeze(MatchingCorrMats(i,j,:)),MatchingPfactorscores); % correlation

    end
    
end

% Sum down columns for composite ROI scores 

% Transpose matrices and add them
BetaMat= BetaMat + BetaMat';
TMat = TMat + TMat';
rMat = rMat + rMat';
% sum down columns
SummedBetas(:,1) = sum(squeeze(BetaMat));
SummedTstats(:,1) = sum(squeeze(TMat));


ABCD = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/ABCD.dtseries.nii');
ABCD.data = zeros(size(ABCD.data,1),1);
ciftilabels = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/Parcels/Parcels_LR.dtseries.nii');

for roi = 1:size(BetaMat,1)
   ABCD.data(ciftilabels.data==roi,1) = SummedTstats(roi,1); 
   %ABCD.data(ciftilabels.data==roi,2) = SummedBetas(roi,2);
   %ABCD.data(ciftilabels.data==roi,3) = SummedBetas(roi,3);
   %ABCD.data(ciftilabels.data==roi,4) = SummedBetas(roi,4);
end
ft_write_cifti_mod('/data/nil-bluearc/GMT/Scott/ABCD/BrainPfactorTstats.dtseries.nii',ABCD)

%% block wise 

close all
% T statistic matrix
%figure; 
imagesc(TMat(NetworksOrdered(:,1),NetworksOrdered(:,1)));colorbar;colormap(jet);caxis([-2 2])
hold on
%draw lines
for n = 1:length(PartitionIdx)
    line([1 length(NetworksOrdered(:,1))],[PartitionIdx(n) PartitionIdx(n)],'Color','k','LineWidth',2.5)
    line([PartitionIdx(n) PartitionIdx(n)],[1 length(NetworksOrdered(:,1))],'Color','k','LineWidth',2.5)
end


