% Merge data tables into one big data table 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/ABCD/Scripts/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/scripts/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/NLAtoolbox/'))

% Read subject list - all the rest subjects we have 
%RestTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/total_subject_list_0327.txt');
RestTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/SubjectLists/total_subjectlist_0508.txt');

RestTable.Properties.VariableNames{1} = 'Subject';
subjects = RestTable.Subject;

 %% Make parcel correlation matrix
% % define input dir 
inputdir = '/data/Daenerys/ABCD/data/abcdbids_output/';
numpools = 7;
%make_abcd_parcel_corrmats(subjects,numpools,inputdir);
% 
 %% Make vertexwise correlation matrix
% inputdir = '/data/Daenerys/ABCD/data/abcdbids_output/';
% outname = '/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/Corrmats_Replication';
% numpools = 7;
% make_abcd_vertex_corrmats(subjects,inputdir,outname);

%% Append correlation matrices & motion information  
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SubjectParcelCorrMats.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFrameTotal.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/Badtmaskidx.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/MissingSubIdx.mat')
for i = 1:size(RestTable,1)
    Thismat = Corrmat(:,:,i);
    Corrmatcell{i,1} = Thismat;
end
RestTable(:,2) = Corrmatcell;
RestTable(:,3) = num2cell(SubjectFD);
RestTable(:,4) = num2cell(SubjectFrameTotal);
RestTable(:,5) = num2cell(MissingSubIdx);
RestTable(:,6) = num2cell(Badtmaskidx);
RestTable.Properties.VariableNames{2} = 'Corrmats';
RestTable.Properties.VariableNames{3} = 'FD';
RestTable.Properties.VariableNames{4} = 'FrameTotal';
RestTable.Properties.VariableNames{5} = 'MissingRestSubjects';
RestTable.Properties.VariableNames{6} = 'Badtmaskidx';

% Scanner info
%ScannerTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/sub_manufac_pairs.txt');
%RestTable = innerjoin(RestTable,ScannerTable,'Keys','Subject');
clear Corrmat
%% Behavior - Discovery set 
%CogTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD_Brain_Cog_Paper/Factorscores_Cognition_4026.txt');

% Discovery CBCL & Total Toolbox Subjects
CBCLTable_g1_discovery = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/G1_Discovery/ABCD_2_0_group1_data_cbcl.txt');
cbcl_total_prob = CBCLTable_g1_discovery.cbcl_scr_syn_totprob_r;
%CBCLTable_g1_discovery(:,2:end) = [];
CBCLTable_g1_discovery = addvars(CBCLTable_g1_discovery,cbcl_total_prob,'After','Subject');


ToolboxTable_g1_discovery = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/G1_Discovery/ABCD_2_0_group1_data_toolboxtotal.txt');
%ToolboxTable_g1_discovery(:,2:7) = [];

% Put them together
% BehavioralTable = innerjoin(CBCLTable_g1_discovery,ToolboxTable_g1_discovery,'Keys','Subject');
% MissingBehavioral = isnan(BehavioralTable.nihtbx_totalcomp_uncorrected); % 0 = not missing 
% BehavioralTable = addvars(BehavioralTable,MissingBehavioral,'After','nihtbx_totalcomp_uncorrected');
MissingBehavioral_g1 = isnan(ToolboxTable_g1_discovery.nihtbx_totalcomp_uncorrected) | isnan(CBCLTable_g1_discovery.cbcl_total_prob); % 0 = not missing 
ToolboxTable_g1_discovery = addvars(ToolboxTable_g1_discovery,MissingBehavioral_g1,'After','nihtbx_totalcomp_uncorrected');
BehavioralTable_g1 = innerjoin(ToolboxTable_g1_discovery,CBCLTable_g1_discovery,'Keys','Subject');
MainTable_g1 = innerjoin(RestTable,BehavioralTable_g1,'Keys','Subject');

%% Behavior - Replication set
% Discovery CBCL & Total Toolbox Subjects
% Discovery CBCL & Total Toolbox Subjects
CBCLTable_g2_replication = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/G2_Replication/ABCD_2_0_group2_data_cbcl.txt');
cbcl_total_prob = CBCLTable_g2_replication.cbcl_scr_syn_totprob_r;
%CBCLTable_g2_replication(:,2:end) = [];
CBCLTable_g2_replication = addvars(CBCLTable_g2_replication,cbcl_total_prob,'After','Subject');


ToolboxTable_g2_replication = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/G2_Replication/ABCD_2_0_group2_data_toolboxtotal.txt');
%ToolboxTable_g1_discovery(:,2:7) = [];

% Put them together
% BehavioralTable = innerjoin(CBCLTable_g1_discovery,ToolboxTable_g1_discovery,'Keys','Subject');
% MissingBehavioral = isnan(BehavioralTable.nihtbx_totalcomp_uncorrected); % 0 = not missing 
% BehavioralTable = addvars(BehavioralTable,MissingBehavioral,'After','nihtbx_totalcomp_uncorrected');
MissingBehavioral_g2 = isnan(ToolboxTable_g2_replication.nihtbx_totalcomp_uncorrected) | isnan(CBCLTable_g2_replication.cbcl_total_prob); % 0 = not missing 
ToolboxTable_g2_replication = addvars(ToolboxTable_g2_replication,MissingBehavioral_g2,'After','nihtbx_totalcomp_uncorrected');
BehavioralTable_g2 = innerjoin(ToolboxTable_g2_replication,CBCLTable_g2_replication,'Keys','Subject');
MainTable_g2 = innerjoin(RestTable,BehavioralTable_g2,'Keys','Subject');

