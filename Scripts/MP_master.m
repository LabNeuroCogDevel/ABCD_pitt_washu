%Master script containing all functions and scripts for RSFC-behavior 
%reproducibility/stability paper 

%% Merge data tables into one big data table 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/ABCD/Scripts/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/scripts/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/NLAtoolbox/'))

% Read subject list - all the rest subjects we have 
RestTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/total_subject_list_0327.txt');
RestTable.Properties.VariableNames{1} = 'Subject';
subjects = RestTable.Subject;

%% Make parcel correlation matrix
% define input dir 
inputdir = '/data/Daenerys/ABCD/data/abcdbids_output/';
numpools = 7;
make_abcd_parcel_corrmats(subjects,inputdir);

%% Make vertexwise correlation matrix
inputdir = '/data/Daenerys/ABCD/data/abcdbids_output/';
outname = '/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/Corrmats_Replication';
numpools = 7;
make_abcd_vertex_corrmats(subjects,inputdir,outname);

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
ScannerTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/sub_manufac_pairs.txt');
RestTable = innerjoin(RestTable,ScannerTable,'Keys','Subject');
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


%% Discovery: FC ~ behavior

% Get variables of interest. 
rmats = cell2mat(MainTable_g1.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g1,1),333);
rmats = permute(rmats,[1,3,2]);

%gen_psych = MainTable_g1.cbcl_scr_syn_totprob_r;
%gen_psych = MainTable.CBCL_F1_1;
gen_cog = MainTable_g1.nihtbx_totalcomp_uncorrected;
%gen_cog = MainTable.ML1;

%FD and threshold by 0.2
FDIdx = MainTable_g1.FD;

%Exclude subjects with missing fMRI rest data
MissingSubs = MainTable_g1.MissingRestSubjects;

% Exclude missing behavioral data
MissingBehavioralData = MainTable_g1.MissingBehavioral_g1;

% Frame total
Frametotal = MainTable_g1.FrameTotal;

% Scanner type 
%scanner = categorical(MainTable_g1.Scanner);
% scanner = dummyvar(scanner);
% scanner = scanner(:,2:3);
scanner = MainTable_g1.Scanner_code;
% symptoms
symptoms = MainTable_g1.cbcl_scr_syn_totprob_r>10;

% Do not include high motion or missing subjects 
rmats_disc = rmats(:,:,FDIdx < 0.2 & Frametotal >= 600 & MissingSubs == 0);% & MissingBehavioralData == 0); %

sublist_disc = MainTable_g1.Subject(FDIdx < 0.2  & Frametotal >= 600 & MissingSubs == 0);

% add whichever factors you want to include here, including interaction term
factors = [gen_cog]; 
% remove high motion subjects and missing rest  and behavioral subs - as above for rmat 
factors = factors(FDIdx < 0.2 & MissingSubs == 0 & MissingBehavioralData == 0 & Frametotal >= 600,:);% & symptoms > 0,:);
scanner = scanner(FDIdx < 0.2 & MissingSubs == 0 & MissingBehavioralData == 0 & Frametotal >= 600);
% zscore
%factors = log(factors);
% Run regression discovery

% Test for interaction? 1 = yes, 0 = no
interaction = 0; 

% Discovery - Run regression for each edge 

[betamat,tmat,pmat,rsq,summed_betas,summed_tstats,summed_pvals,corrs] = run_regression(rmats_disc,factors,interaction);

%% Replication: FC ~ behavior

% Get variables of interest. 
rmats = cell2mat(MainTable_g2.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g2,1),333);
rmats = permute(rmats,[1,3,2]);

%gen_psych = MainTable_g2.cbcl_scr_syn_totprob_r;
%gen_psych = MainTable.CBCL_F1_1;
gen_cog = MainTable_g2.nihtbx_totalcomp_uncorrected;
%gen_cog = MainTable.ML1;

%FD and threshold by 0.2
FDIdx = MainTable_g2.FD;

%Exclude subjects with missing fMRI rest data
MissingSubs = MainTable_g2.MissingRestSubjects;

% Exclude missing behavioral data
MissingBehavioralData = MainTable_g2.MissingBehavioral_g2;

% Frame total
Frametotal = MainTable_g2.FrameTotal;

% Scanner type 
% scanner = categorical(MainTable_g2.Scanner);
% scanner = dummyvar(scanner);
% scanner = scanner(:,2:3);
scanner = MainTable_g2.Scanner_code;

% symptoms
symptoms = MainTable_g2.cbcl_scr_syn_totprob_r>10;

% Do not include high motion or missing subjects 
rmats_rep = rmats(:,:,FDIdx < 0.2  & Frametotal >= 600 & MissingSubs == 0);% & MissingBehavioralData == 0); %
sublist_rep = MainTable_g2.Subject(FDIdx < 0.2  & Frametotal >= 600);
%zscore total scores
%gen_psych = zscore(gen_psych);
%gen_cog = zscore(gen_cog);
%sem = MainTable.SES_composite;

% add whichever factors you want to include here, including interaction term
factors = [gen_cog]; 
% remove high motion subjects and missing rest  and behavioral subs - as above for rmat 
factors = factors(FDIdx < 0.2 & MissingSubs == 0 & MissingBehavioralData == 0 & Frametotal >= 600,:);% & symptoms > 0,:);

% zscore
%factors = log(factors);
% Replication - run regression for each edge

[betamat,tmat,pmat,rsq,summed_betas,summed_tstats,summed_pvals,corrs] = run_regression(rmats_rep,factors,interaction);

%% Plotting 
% Reorder ROIs from gordon laumann parcellation & create partition index
reorder_gordon_laumann_parcels
% Get ROI order to pass into plotting function 
roi_order = NetworksOrdered(:,1);
%Plot it
plot_adj_matrix(STD_rmatsrep_ge,roi_order,partitionidx);
%tmean=mean(nanmean(tmat));
%mini=tmean - 2*(mean(std(tmat)));
%maxi=tmean + 2*(mean(std(tmat)));
%caxis([mini maxi]); 
caxis([.1 .25]);
colormap(hot);

%% Blockwise iterative RSFC-behavior split-half reliability 

% Obtain rmats and factor from above for input to function 
allmats = cat(3,rmats_disc,rmats_rep);
allmats_ordered = allmats(roi_order,roi_order,:);
factor = [factors;factors_rep];

[NetworkCorrs,NetworkPvals] = abcd_blockwise_correlation_iterative_reliability(allmats_ordered,factor,partitionidx,25,1000,1,7);

% for q = 1:size(NetworkPvals,1)
% output(q,1) = 100*(length(find(NetworkPvals(q,:)>.05)))/size(NetworkPvals,2);
% end
% type1errprobability = sort(output,'descend');
% 
% % Plot it 
% plot(type1errprobability,'k','LineWidth',3)
% xticks([4:8:85])
% xticklabels({'100','300','500','700','900','1100','1300','1500','1700','1900','2100'})
% ylabel('Probability of Committing Type II Error')
% xlabel('Number of Subjects')
% set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
% title('RSFC NIH Toolbox DMN/DAN Brain Behavior Correlation')

% Get +/- 2 SDs and plot using shaded error bars 
figure;
hold on
%CBCL
Mean_netcorr = mean(NetworkCorrs_cbcl,2)';
lineProps.col = {[1 127/255 80/255]};
stdcorr = std(NetworkCorrs_cbcl').*2;
h = mseb([1:1:85],Mean_netcorr,stdcorr,lineProps,1);

%Tlbx
Mean_netcorr = mean(NetworkCorrs_cog,2)';
lineProps.col = {[32/255 178/255 170/255]};
stdcorr = std(NetworkCorrs_cog').*2;
k = mseb([1:1:85],Mean_netcorr,stdcorr,lineProps,1);


set(gca,'FontWeight','bold','FontSize',14)
xlim([1 85]); ylim([-.5 .5])
xticks([0 12:12:85]) % every 300 subjects; each bin = 25 subs
xticklabels({'25'; '300';'600';'900';'1200';'1500';'1800';'2100'});
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Brain Behavior Correlation','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)