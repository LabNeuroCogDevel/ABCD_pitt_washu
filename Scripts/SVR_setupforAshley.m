%% ABCD Brain Behavior 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts'))
% Set input variables. 
%Sub list
subjectlist = '/data/nil-bluearc/GMT/ABCD/Initial_Subjectlist_287.txt';
% Base directory to cifti time series 
InputDir = '/data/nil-bluearc/GMT/ABCD/initial_ABCD_data/';
% Path to factor scores - assumes format SID ML1 ML2 ML3 ML4 SID
FactorScoreFile = '/data/nil-bluearc/GMT/Scott/ABCD/factorscores_exploratory267.txt';

%% Made group averaged dconn
% Read in gordon parcels 
ciftilabels = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/Parcels/Parcels_LR.dtseries.nii');
% Read in sub list
subjects = textread(subjectlist,'%s','delimiter','\n');

for s = 1:size(subjects,1)
    disp(['Subject ' num2str(s)])
    AverageTS = [];
    if s ~= 4
        subject = subjects{s};
        a = dir([InputDir subject '/MNINonLinear/Results/*rfMRI_REST*']);
        b = dir([InputDir subject '/MNINonLinear/Results/*rfMRI_REST_*']);
        n = numel(a) - numel(b);
        FD = [];
        for run = 1:n
            ThisFD = load([InputDir subject '/MNINonLinear/Results/rfMRI_REST' num2str(run) '/FNL_preproc_v2/FD.mat']);
            data = zeros(1+length(ThisFD.FD),1);
            data(1) = 0; data(2:end) = ThisFD.FD;
            FD = [FD; data];
        end
        SubjectFD(s,1) = mean(FD(2:end));
        FD(FD>=.2)=0; FD(FD<.2 & FD > 0) = 1; ThisFD = sum(FD);
        SubjectFrameTotal(s,1) = ThisFD;
    else
        SubjectFD(s,1) = 1;
        SubjectFrameTotal(s,1) = 0;
    end
end

load('/data/nil-bluearc/GMT/Scott/ABCD_DCN_CognitionPaper/SubjectCorrMats.mat','Corrmat')

%% Read in behavioral data and align subjects between brain & behavior
[CogSublist, ML_one, ML_two, ML_four, ML_three, Subkey] = textread(FactorScoreFile,'%s%s%s%s%s%s');
%CogSublist = textread('/data/nil-bluearc/GMT/Scott/ABCD/factorscores_exploratory267.txt','%s','delimiter','\n');
CogSublist = CogSublist(2:end,1);
ML_one = ML_one(2:end,1);
ML_two = ML_two(2:end,1);
ML_three = ML_three(2:end,1);
ML_four = ML_four(2:end,1);

for n = 1:size(ML_one,1)
    Factors(n,1) = str2num(ML_one{n,1});
    Factors(n,2) = str2num(ML_two{n,1});
    Factors(n,3) = str2num(ML_three{n,1});
    Factors(n,4) = str2num(ML_four{n,1});
end

BrainSublist = textread(subjectlist,'%s','delimiter','\n');

for s = 1:size(BrainSublist,1)
    ThisSub = BrainSublist{s}; 
    BrainSublist{s,1} = ThisSub(5:end);
end

% Line up Behavioral subjects with Brain Subjects 
MatchingIdx = zeros(size(BrainSublist,1),1);
for sub = 1:size(CogSublist,1)
    ThisCogSub = CogSublist{sub};
    for BrainSub = 1:size(BrainSublist,1)
        ThisBrainSub = BrainSublist{BrainSub};
        if all(ThisCogSub == ThisBrainSub)
            MatchingIdx(BrainSub,1) = 1;
        end
    end
end
MatchingBrainSublist = BrainSublist(MatchingIdx==1);
ComboList = [];
ComboList = [CogSublist MatchingBrainSublist];

% Check to ensure they are lined up 
List = zeros(size(CogSublist,1),1);
for s = 1:size(ComboList,1)
   ThisCogSub = CogSublist{s};
   ThisBrainSub = MatchingBrainSublist{s};
   if ThisCogSub == ThisBrainSub
       List(s,1) = 1;
   end
end
sum(List)

if sum(List) ~= size(ComboList,1)
    disp('Mismatch is brain/behavior subjust alignment!!')
end
MatchingSubsFD=SubjectFD(MatchingIdx==1);
FDIdx=find(MatchingSubsFD<=.2);
% Remove Subjects from brain data not in behavioral data
MatchingCorrMats = Corrmat(:,:,MatchingIdx==1);
MatchingCorrMats=MatchingCorrMats(:,:,FDIdx);
Factors = Factors(FDIdx,:);

%% Get behavioral measures (NIH tlbx)
NIHToolbox=readtable('/data/nil-bluearc/GMT/Scott/ABCD/ABCD_NIHtoolboxScores.csv');


ABCD_Subjects_All = textread('/data/nil-bluearc/GMT/Scott/ABCD/ABCD_AllSubjects.txt','%s');
ABCD_Subjects_Partial = textread('/data/nil-bluearc/GMT/ABCD/Initial_Subjectlist_287.txt','%s');

% Remove _ in big list IDs
for s = 1:size(ABCD_Subjects_All)
    ThisSub = ABCD_Subjects_All{s};
    SubListBig(s,:) = ThisSub(1,[1:4 6:end]);
end
SubjectMatchIdx = zeros(size(SubListBig,1),1);
% Get index of each subject in big list 
for s = 1:size(ABCD_Subjects_Partial,1)
    subject = ABCD_Subjects_Partial{s};
    subject = subject(5:end);
    for t = 1:size(ABCD_Subjects_All)
        if all(SubListBig(t,:) == subject)
            SubjectMatchIdx(t,1) = 1;
        end
    end
end
SubjectSubset = ABCD_Subjects_All(SubjectMatchIdx==1);

for vars = 3:size(NIHToolbox,2)
    Input = table2array(NIHToolbox(:,3));
    ToolboxTotal=Input(SubjectMatchIdx==1);
    
end

ToolboxTotalSubset = ToolboxTotal(FDIdx,1);
ToolboxTotalMatchingIdx = [];
for n = 1:size(ToolboxTotalSubset,1)
    Thisdata = ToolboxTotalSubset{n};
    if Thisdata(1)=='#'
        ToolboxTotalMatchingIdx(n,1) = 0;
    else
        ToolboxTotalMatchingIdx(n,1) = 1;
    end
end

NewMatchingCorrmats = MatchingCorrMats(:,:,ToolboxTotalMatchingIdx==1);
MatchingToolboxTotals = ToolboxTotalSubset(ToolboxTotalMatchingIdx==1);
for n = 1:size(MatchingToolboxTotals,1)
   NewMatchingToolboxTotals(n,1) = str2num(MatchingToolboxTotals{n});
end

save('/data/nil-bluearc/GMT/Scott/ABCD/Corrmats_Ashley','NewMatchingCorrmats')
save('/data/nil-bluearc/GMT/Scott/ABCD/ToolboxTotal_Ashley','NewMatchingToolboxTotals')


