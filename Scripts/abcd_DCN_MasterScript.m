% Merge data tables into one big data table 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/ABCD/Scripts/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/scripts/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/NLAtoolbox/'))

% Read subject list - all the rest subjects we have 
RestTable_DCN = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/total_subject_list_0327.txt');
RestTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/SubjectLists/total_subjectlist_0508.txt');

RestTable.Properties.VariableNames{1} = 'Subject';
RestTable_DCN.Properties.VariableNames{1} = 'Subject';
subjects = RestTable.Subject;

%% Make parcel correlation matrix
% define input dir 
inputdir = '/data/Daenerys/ABCD/data/abcdbids_output/';
numpools = 7;
make_abcd_parcel_corrmats(subjects,numpools,inputdir);

%% Make vertexwise correlation matrix
inputdir = '/data/Daenerys/ABCD/data/abcdbids_output/';
outname = '/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/Corrmats_Replication';
numpools = 7;
make_abcd_vertex_corrmats(subjects,inputdir,outname);

%% Append correlation matrices & motion information  
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SubjectParcelCorrMats.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD_Incframes.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFrameTotal.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/Badtmaskidx.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/MissingSubIdx.mat')
for i = 1:size(RestTable,1)
    Thismat = Corrmat(:,:,i);
    Corrmatcell{i,1} = Thismat;
end
RestTable(:,2) = Corrmatcell;
RestTable(:,3) = num2cell(SubjectFD);
RestTable(:,4) = num2cell(SubjectFD_Incframes);
RestTable(:,5) = num2cell(SubjectFrameTotal);
RestTable(:,6) = num2cell(MissingSubIdx);
RestTable(:,7) = num2cell(Badtmaskidx);
RestTable.Properties.VariableNames{2} = 'Corrmats';
RestTable.Properties.VariableNames{3} = 'FD';
RestTable.Properties.VariableNames{4} = 'FD_Incframes';
RestTable.Properties.VariableNames{5} = 'FrameTotal';
RestTable.Properties.VariableNames{6} = 'MissingRestSubjects';
RestTable.Properties.VariableNames{7} = 'Badtmaskidx';

% Scanner info
%ScannerTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/sub_manufac_pairs.txt');
%RestTable = innerjoin(RestTable,ScannerTable,'Keys','Subject');
RestTable_DCN = innerjoin(RestTable_DCN,RestTable,'Keys','Subject');
RestTable = RestTable_DCN;
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

gen_psych = MainTable_g1.cbcl_scr_syn_totprob_r;
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
%scanner = MainTable_g1.Scanner_code;
% symptoms
symptoms = MainTable_g1.cbcl_scr_syn_totprob_r>10;

% Do not include high motion or missing subjects 
rmats_disc = rmats(:,:, MissingSubs == 0); %

sublist_disc = MainTable_g1.Subject(FDIdx < 0.2  & Frametotal >= 600 & MissingSubs == 0);

% add whichever factors you want to include here, including interaction term
factors = [gen_cog];% FDIdx]; 
% remove high motion subjects and missing rest  and behavioral subs - as above for rmat 
factors = factors(FDIdx < 0.2 & MissingSubs == 0 & MissingBehavioralData == 0 & Frametotal >= 600,:);% & symptoms > 0,:);
%scanner = scanner(FDIdx < 0.2 & MissingSubs == 0 & MissingBehavioralData == 0 & Frametotal >= 600);
% zscore
%factors = log(factors);
% Run regression discovery

% Test for interaction? 1 = yes, 0 = no
interaction = 0; 

% Discovery - Run regression for each edge 

[betamat_disc,tmat_disc,pmat_disc,rsq,summed_betas_disc,summed_tstats_disc,summed_pvals_disc,corrs_disc] = run_regression(rmats_disc,factors,interaction);

%% Scanner regession 
[rsqdiff_disc,mean_rfiddsq_disc] = run_regression_scanner(rmats_disc,factors,scanner);

%% Replication: FC ~ behavior

% Get variables of interest. 
rmats = cell2mat(MainTable_g2.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g2,1),333);
rmats = permute(rmats,[1,3,2]);

gen_psych = MainTable_g2.cbcl_scr_syn_totprob_r;
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
%scanner = MainTable_g2.Scanner_code;

% symptoms
symptoms = MainTable_g2.cbcl_scr_syn_totprob_r>10;

% Do not include high motion or missing subjects 
rmats_rep = rmats(:,:,FDIdx < 0.2  & Frametotal >= 600 & MissingSubs == 0 & MissingBehavioralData == 0); %
sublist_rep = MainTable_g2.Subject(FDIdx < 0.2  & Frametotal >= 600);
%zscore total scores
%gen_psych = zscore(gen_psych);
%gen_cog = zscore(gen_cog);
%sem = MainTable.SES_composite;

% add whichever factors you want to include here, including interaction term
factors_rep = [gen_cog FDIdx]; 
% remove high motion subjects and missing rest  and behavioral subs - as above for rmat 
factors_rep = factors_rep(FDIdx < 0.2 & MissingSubs == 0 & MissingBehavioralData == 0 & Frametotal >= 600,:);% & symptoms > 0,:);

% zscore
%factors = log(factors);
% Replication - run regression for each edge

[betamat_rep,tmat_rep,pmat_rep,rsq,summed_betas_rep,summed_tstats_rep,summed_pvals_rep,corrs_rep] = run_regression(rmats_rep,factors_rep,interaction);

%% Plotting 
% Reorder ROIs from gordon laumann parcellation & create partition index
reorder_gordon_laumann_parcels
% Get ROI order to pass into plotting function 
roi_order = NetworksOrdered(:,1);
%Plot it
plot_adj_matrix(tmat_rep(:,:,1),roi_order,partitionidx);

caxis([-.4 1]);
colormap(jet);
%% Save
outname = 'CBCL_pfactor_t.mat';
out = tmat;
save(['/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/' outname],'out');
outname = 'CBCL_TLBX_SEM_tlbxloading_t.mat';
out = tmat(:,:,2);
save(['/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/' outname],'out');
outname = 'CBCL_TLBX_SEM_semloading_t.mat';
out = tmat(:,:,3);
save(['/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/' outname],'out');
%% Write file 
ABCD = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/ABCD.dtseries.nii');
ABCD.data = zeros(size(ABCD.data,1),3);
ciftilabels = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/Parcels/Parcels_LR.dtseries.nii');
outprefix = 'Toolboxtotal_Disc_Rep';
for roi = 1:size(rmats,1)
    ABCD.data(ciftilabels.data==roi,1) = summed_tstats_disc(roi,1);
    ABCD.data(ciftilabels.data==roi,2) = summed_tstats_rep(roi,1); 
    ABCD.data(ciftilabels.data==roi,3) = summed_tstats_disc(roi,1)-summed_tstats_rep(roi,1);
end
ABCD.time = [1:3];
cifti = ['/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/' outprefix '.dtseries.nii'];
ft_write_cifti_mod(cifti,ABCD)
%cd('/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/')
%wb_custom_colorbar(cifti,[-100 100])

%% Enrichment
% set params
IM.name = 'CBCL_Discovery';
make_CLUT_power17;
CLUT = ColorLUT([1:3 5 7:12 15 16],:);
IM.cMap = CLUT;
IM.Nets = {'DMN';'Vis';'FPN';'DAN';'VAN';'SAL';'CON';'SMd';'SMl';'Aud';'ParMem';'CTXT'};
read_parcelinfo; % generates a table to add Gordon ROI info to the IM struct 
CrapNetwork = ParcelInfo.NetworkID==13;
ParcelInfo(CrapNetwork,:)=[];
IM.ROIxyz = [ParcelInfo.X ParcelInfo.Y ParcelInfo.Z];
IM.key = [ParcelInfo.ParcelID ParcelInfo.NetworkID];
% Reorder ROIs from gordon laumann parcellation & create partition index
reorder_gordon_laumann_parcels
NetworksOrdered(CrapNetwork,:) = [];
IM.order = NetworksOrdered(:,1);
fc = rmats_disc;
fc(CrapNetwork,:,:)=[];
fc(:,CrapNetwork,:)=[];
% Define a behavior name 
Bxname='CBCL_disc';
Bx = factors;
%% Run pipeline
Enrichment_pipeline

%% STD and CoV
MeanCorrmat_disc = mean(rmats_disc,3);
STDCorrmat_disc = std(rmats_disc,0,3);
CoeffofVar_disc = STDCorrmat_disc ./ MeanCorrmat_disc;
MeanCorrmat_rep = mean(rmats_rep,3);
STDCorrmat_rep = std(rmats_rep,0,3);
CoeffofVar_rep = STDCorrmat_rep ./ MeanCorrmat_rep;

% By brain region 
Region_std_disc = mean(STDCorrmat_disc)';
Region_std_rep = mean(STDCorrmat_rep)';
corr(Region_std_disc,Region_std_rep)
% Write file 
ABCD = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/ABCD.dtseries.nii');
ABCD.data = zeros(size(ABCD.data,1),3);
ciftilabels = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/Parcels/Parcels_LR.dtseries.nii');
outprefix = 'STD_Disc_Rep';
for roi = 1:size(rmats_disc,1)
    ABCD.data(ciftilabels.data==roi,1) = Region_std_disc(roi,1);
    ABCD.data(ciftilabels.data==roi,2) = Region_std_rep(roi,1); 
    ABCD.data(ciftilabels.data==roi,3) = Region_std_disc(roi,1)-Region_std_rep(roi,1);
end
ABCD.time = [1:3];
cifti = ['/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/' outprefix '.dtseries.nii'];
ft_write_cifti_mod(cifti,ABCD)


%% correlation between discovery and replication brain/behavior 
%edgewise
tmat_disc_c = tmat_disc_cogfdmain;
tmat_rep_c = tmat_disc_cogonly;
uidx = find(triu(tmat_disc_c,1));
Disc_Rep_corr_edgewise = FisherTransform(corr(tmat_disc_c(uidx),tmat_rep_c(uidx)););

% Blockwise (by network)
MeanMatrix_disc = abcd_quantify_blocks(tmat_disc_c,partitionidx,roi_order);
MeanMatrix_rep = abcd_quantify_blocks(tmat_rep_c,partitionidx,roi_order);
MeanMatrix_disc = MeanMatrix_disc(1:12,1:12);
MeanMatrix_rep = MeanMatrix_rep(1:12,1:12);
uidx = find(triu(MeanMatrix_disc));
Disc_Rep_corr_block = FisherTransform(corr(MeanMatrix_disc(uidx),MeanMatrix_rep(uidx)));
%% Split half reliability BrBx

% Get 1st half of data Brain behavior corrs (FHBrBx)
FHBrBx = zeros(size(rmats_rep,1));
for i = 1:(size(rmats_rep,1)-1)
    disp(['On ROI ' num2str(i)])
    for j = i+1:size(rmats_rep,1)
        FHBrBx(i,j)  = corr(squeeze(rmats_rep(i,j,:)),factors);
    end
end
FHBrBx = FHBrBx + FHBrBx';
corrmats = rmats_disc;
%gen_cog = MainTable_g1.nihtbx_totalcomp_uncorrected;
gen_psych = MainTable_g1.cbcl_scr_syn_totprob_r;
FDIdx = MainTable_g1.FD;
MissingSubs = MainTable_g1.MissingRestSubjects;
MissingBehavioralData = MainTable_g1.MissingBehavioral_g1;
Frametotal = MainTable_g1.FrameTotal;
behavior = gen_psych; 
behavior = behavior(FDIdx < 0.2 & MissingSubs == 0 & MissingBehavioralData == 0 & Frametotal >= 600,:);% & symptoms > 0,:);


[ReliabilityCoeffs_edgewise, ReliabilityCoeffs_network, ReliabilityCoeffs_blocks] = abcd_xsub_brainbehav_reliability(FHBrBx,corrmats,behavior,25,100,partitionidx,roi_order);

%% Plot
Thiscurve_tlbx = squeeze(ReliabilityCoeffs_blocks(10,10,:,:));
Thiscurve_cbcl = squeeze(ReliabilityCoeffs_blocks_CBCL(1,9,:,:));
figure;
hold on
%CBCL
Mean_edge_reliability = mean(Thiscurve_cbcl,2);
lineProps.col = {[1 127/255 80/255]};
Semedge_reliability = std(Thiscurve_cbcl').*2;
h = mseb([1:1:44],Mean_edge_reliability+.04,Semedge_reliability,lineProps,1);

%Tlbx
Mean_edge_reliability = mean(Thiscurve_tlbx,2);
lineProps.col = {[32/255 178/255 170/255]};
Semedge_reliability = std(Thiscurve_tlbx').*2;
k = mseb([1:1:44],Mean_edge_reliability+.04,Semedge_reliability,lineProps,1);


set(gca,'FontWeight','bold','FontSize',14)
xlim([1 44])
xticklabels({'250','500','750','1000'});
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Split-half Pearson r','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)



%% Average Dconn & Infomap: discovery 

% Load Distance matrix
dmat = '/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/MSC01/fsaverage_LR32k/MSC01_distmat_surf_geodesic_vol_euclidean_MNInlwarp_xhem1000_uint8.mat';
% Load Dconn
dconn = '/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/Corrmats_Replication.dconn.nii';
% params
xdist = 30;
thresholdarray = [.001:.001:.01 .02:.01:.05];
%thresholdarray = [.001:.001:.007];
makebinary = 0;
numpools = 7;
structure_indices = [];
outdir = '/data/nil-bluearc/GMT/Scott/ABCD/Infomap/Replication';

% Infomap
Run_Infomap_2015_SM(dconn,dmat,xdist,thresholdarray,makebinary,outdir,numpools,structure_indices)
clear dconn dmat %dmat
% Clean up 
simplified = modify_clrfile('simplify',[outdir '/rawassn.txt'],400);
regularized = regularize(simplified);
cifti_struct = '/data/nil-bluearc/GMT/Evan/MSC/Analysis_V2/Networks_template.dscalar.nii';
cifti_struct = ft_read_cifti_mod(cifti_struct);
Datasize = size(cifti_struct.data,1);
cifti_struct.data = regularized;
cifti_struct.data(59413:Datasize,:) = 0;
cifti_struct.dimord = 'pos_time';
cifti_struct.time=[1:1:length(thresholdarray)];
ft_write_cifti_mod([outdir '/GroupAvg_rawassn_minsize400_regularized'],cifti_struct);


%% MDS plots 
Var_discovery = MainTable_g1.demo_gender_id_v2_l_ToolboxTable_g1_discovery;
Var_discovery = Var_discovery(MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0);
rmats = cell2mat(MainTable_g1.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g1,1),333);
rmats = permute(rmats,[1,3,2]);
rmats_disc = rmats(:,:,MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0);

Var_replication = MainTable_g2.demo_gender_id_v2_l_ToolboxTable_g2_replication;
Var_replication = Var_replication(MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0);
rmats = cell2mat(MainTable_g2.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g2,1),333);
rmats = permute(rmats,[1,3,2]);
rmats_rep = rmats(:,:,MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0);

[Y_disc E_disc h_disc] = cmdscale_mat_sm(rmats_disc,Var_discovery,'euclidean',groupnames);
[Y_rep E_rep h_rep] = cmdscale_mat_sm(rmats_rep,Var_replication,'euclidean',groupnames);

% Quantify via PCA 
[coeff, score, latent, tsqaured, explain, mu] = pca(Y_disc);
First_ten_variance_disc = sum(explain(1:10));
[coeff, score, latent, tsqaured, explain, mu] = pca(Y_rep);
First_ten_variance_rep = sum(explain(1:10));

% Site similarity matrix 
[similarity_matrix, partitions, order] = abcd_site_similarity(rmats_disc,Var_discovery);

% Plotting
plot_adj_matrix(similarity_matrix,order,partitions);
caxis([0 1]);
colormap(parula);

% Quantify by block
MeanMatrix = abcd_quantify_blocks(similarity_matrix,partitions,order);

% Get stats on block
[t_disc,p_disc,d_disc] = abcd_block_stats(similarity_matrix,partitions(1),order);
[t_rep,p_rep,d_rep] = abcd_block_stats(similarity_matrix,partitions(1),order);

% MDS of two most dissimilar sites
matrix = MeanMatrix;
minval = min(min(matrix));
[row, col] = find(matrix==minval);
row = row(1); col = col(1);
subidx = order(partitions(row-1)+1:partitions(row)); 
subidx= [subidx;order(partitions(col-1)+1:partitions(col))];
group = ones(length(subidx),1);
g1length = length(order(partitions(row-1)+1:partitions(row))); 
group(g1length+1:end) = group(g1length+1:end)+1;
[Y E h] = cmdscale_mat_sm(rmats_rep(:,:,subidx),group,'euclidean');

% centroid circle
Thisvar = Var_replication;
figure; hold on
for i = 1:numel(unique(Thisvar))
    xy_centroid = mean([Y(Thisvar==i,1),Y(Thisvar==i,2)]);
    scatter(xy_centroid(1),xy_centroid(2),100,colors(i,:),'filled')
    distances_mds = .5*mean(pdist([Y(Thisvar==i,1),Y(Thisvar==i,2)]));
    viscircles(xy_centroid,distances_mds,'Color',colors(i,:),'LineWidth',3);
end
xlabel('MDS 1','FontWeight','bold','FontSize',14)
ylabel('MDS 2','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
xlim([-15 15])
ylim([-15 15])
legend(groupnames,'FontWeight','bold','FontSize',14)

%% Scanner effects by distance
read_gordon_coordinates
xyz = table2array(ParcelsNetworkIds);
parceldmat = pdist(xyz);
parceldmat = squareform(parceldmat);
Thismat = corrs_rep;
uidx = find(triu(Thismat,1));
figure;
plot(parceldmat(uidx),Thismat(uidx),'.','Color',[0 0 0],'MarkerSize',10);
set(gca,'FontWeight','bold','FontSize',14)
xlabel('Euclidean Distance','FontWeight','bold','FontSize',14)
ylabel('RSFC/Scanner F-statistic','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
xlim([0 180])
%ylim([-.4 .4])

%% Scanner effects ANOVA
% Discovery set 
rmats = cell2mat(MainTable_g1.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g1,1),333);
rmats = permute(rmats,[1,3,2]);
rmats_disc = rmats(:,:,MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0);

groups = MainTable_g1.Scanner;
groups = groups(MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0);
Fstats = zeros(size(rmats_disc,1));
pv = zeros(size(rmats_disc,1));
for i = 1:size(rmats_disc,1)
    for j = i+1:size(rmats_disc,2)
        [pv(i,j) tbl] = anova1(squeeze(rmats_disc(i,j,:)),groups,'off');
        Fstats(i,j) = tbl{2,5};
    end
end
Fstats_disc = Fstats + Fstats';
pv = pv + pv';

% Quantify significant results 
uidx = find(triu(pv,1));
upperpv = pv(uidx);
sigedges_disc = length(upperpv(upperpv<0.001)) / length(upperpv);

% Replication
rmats = cell2mat(MainTable_g2.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g2,1),333);
rmats = permute(rmats,[1,3,2]);
rmats_rep = rmats(:,:,MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0);

groups = MainTable_g2.Scanner;
groups = groups(MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0);
Fstats = zeros(size(rmats_rep,1));
pv = zeros(size(rmats_rep,1));
for i = 1:size(rmats_rep,1)
    for j = i+1:size(rmats_rep,2)
        [pv(i,j) tbl] = anova1(squeeze(rmats_rep(i,j,:)),groups,'off');
        Fstats(i,j) = tbl{2,5};
    end
end
Fstats_rep = Fstats + Fstats';
pv = pv + pv';

Fstats_sum_disc = mean(Fstats_disc,2);
Fstats_sum_rep = mean(Fstats_rep,2);

% Quantify significant results 
uidx = find(triu(pv,1));
upperpv = pv(uidx);
sigedges_rep = length(upperpv(upperpv<0.001)) / length(upperpv);

% Write file 
ABCD = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/ABCD.dtseries.nii');
ABCD.data = zeros(size(ABCD.data,1),2);
ciftilabels = ft_read_cifti_mod('/data/nil-bluearc/GMT/Scott/Parcels/Parcels_LR.dtseries.nii');
outprefix = 'Fstats_scanner_mean_Disc_Rep';
for roi = 1:size(rmats_disc,1)
    ABCD.data(ciftilabels.data==roi,1) = Fstats_sum_disc(roi,1);
    ABCD.data(ciftilabels.data==roi,2) = Fstats_sum_rep(roi,1); 
end
ABCD.time = [1:2];
cifti = ['/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/' outprefix '.dtseries.nii'];
ft_write_cifti_mod(cifti,ABCD)


%% Scanner -> GE/Philips collapsed; siemens separate 

rmats = cell2mat(MainTable_g1.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g1,1),333);
rmats = permute(rmats,[1,3,2]);
rmats_disc = rmats(:,:,MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0);

groups = MainTable_g1.Scanner_code;
groups = groups(MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0);
groups(groups==3) = 2;

% Replication
rmats = cell2mat(MainTable_g2.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g2,1),333);
rmats = permute(rmats,[1,3,2]);
rmats_rep = rmats(:,:,MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0);

groups = MainTable_g2.Scanner_code;
groups = groups(MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0);
groups(groups==3) = 2;


corrs_rep = zeros(size(rmats_disc,1),size(rmats_disc,2));
for i = 1:size(rmats_disc,1)
    for j = i+1:size(rmats_disc,2)
        corrs_rep(i,j) = corr(squeeze(rmats_rep(i,j,:)),groups);
    end
end
corrs_rep = corrs_rep + corrs_rep';




%% PCA
%Discovery 
PCA_disc = load('/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/discovery_pcaCoeff_roiXroiXcomp_All.mat');

%Replication
PCA_rep = load('/data/nil-bluearc/GMT/Scott/ABCD/DCNpaper/replication_pcaCoeff_roiXroiXcomp_All.mat');


% PCA reliability (disc/rep comparison)
for c = 1:10
    Component = PCA_disc.roi_comps(:,:,c); 
    Component(isnan(Component))=0;
    Component_disc(:,:,c) = Component + Component';
    
    uidx = find(triu(squeeze(Component_disc(:,:,c)),1));
    
    Component = PCA_rep.roi_comps(:,:,c); 
    Component(isnan(Component))=0;
    Component_rep(:,:,c) = Component + Component';
    
    ThisD = squeeze(Component_disc(:,:,c));
    ThisR = squeeze(Component_rep(:,:,c));
    
    Component_reliability(c,1) = corr(ThisD(uidx),ThisR(uidx));
    if Component_reliability(c,1) < 0 % flip sign to match discovery set
        Component_rep(:,:,c) = Component_rep(:,:,c).*-1;
    end
end

reorder_gordon_laumann_parcels
% Get ROI order to pass into plotting function 
roi_order = NetworksOrdered(:,1);
%Plot it
for c = 1:10
    %figure;
    
    plot_adj_matrix(Component_rep(:,:,c),roi_order,partitionidx);
    colormap(jet)
    caxis([-.02 .02])
    title(['Component ' num2str(c)])
    set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
    xticks([])
    yticks([])
    
end

%% Default mode / Dorsal Attn Correlations in subsets of subjects 

% Obtain rmats and factor from above for input to function 
allmats = cat(3,rmats_disc,rmats_rep);
allmats_ordered = allmats(roi_order,roi_order,:);
factors = [factors;factors_rep];

[NetworkCorrs,NetworkPvals] = abcd_blockwise_correlation_iterative_reliability(allmats_ordered,factors,partitionidx,25,1000,1,1);

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

%% Supplement 
rmats = cell2mat(MainTable_g1.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g1,1),333);
rmats = permute(rmats,[1,3,2]);
rmats_disc = rmats(:,:,MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0);

groups = MainTable_g1.Scanner_code;
groups = groups(MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0);

rmats_disc_siemens = rmats_disc(:,:,groups==1);
rmats_disc_philips = rmats_disc(:,:,groups==2);
rmats_disc_ge = rmats_disc(:,:,groups==3);

STD_rmatsdisc_siemens = std(rmats_disc_siemens,0,3);
STD_rmatsdisc_philips = std(rmats_disc_philips,0,3);
STD_rmatsdisc_ge = std(rmats_disc_ge,0,3);

%[p tbl] = anova1([STD_rmatsdisc_siemens(uidx),STD_rmatsdisc_philips(uidx),STD_rmatsdisc_ge(uidx)]);

% effect size 
M1 = mean(STD_rmatsdisc_siemens(uidx));
M2 = mean(STD_rmatsdisc_philips(uidx));
M3 = mean(STD_rmatsdisc_ge(uidx));
SD1 = std(STD_rmatsdisc_siemens(uidx));
SD2 = std(STD_rmatsdisc_philips(uidx));
SD3 = std(STD_rmatsdisc_ge(uidx));

% Siemens / Philips
sdpooled = sqrt((SD1^2 + SD2^2) / 2);
d_Siem_Phil_disc = (M1-M2) ./ sdpooled;

% Siemens / Philips
sdpooled = sqrt((SD1^2 + SD3^2) / 2);
d_Siem_GE_disc = (M1-M3) ./ sdpooled;

% Siemens / Philips
sdpooled = sqrt((SD2^2 + SD3^2) / 2);
d_Phil_GE_disc = (M2-M3) ./ sdpooled;


% Replication
rmats = cell2mat(MainTable_g2.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g2,1),333);
rmats = permute(rmats,[1,3,2]);
rmats_rep = rmats(:,:,MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0);

groups = MainTable_g2.Scanner_code;
groups = groups(MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0);

rmats_rep_siemens = rmats_rep(:,:,groups==1);
rmats_rep_philips = rmats_rep(:,:,groups==2);
rmats_rep_ge = rmats_rep(:,:,groups==3);

STD_rmatsrep_siemens = std(rmats_rep_siemens,0,3);
STD_rmatsrep_philips = std(rmats_rep_philips,0,3);
STD_rmatsrep_ge = std(rmats_rep_ge,0,3);

%[p tbl] = anova1([STD_rmatsdisc_siemens(uidx),STD_rmatsdisc_philips(uidx),STD_rmatsdisc_ge(uidx)]);

% effect size 
M1 = mean(STD_rmatsrep_siemens(uidx));
M2 = mean(STD_rmatsrep_philips(uidx));
M3 = mean(STD_rmatsrep_ge(uidx));
SD1 = std(STD_rmatsrep_siemens(uidx));
SD2 = std(STD_rmatsrep_philips(uidx));
SD3 = std(STD_rmatsrep_ge(uidx));

% Siemens / Philips
sdpooled = sqrt((SD1^2 + SD2^2) / 2);
d_Siem_Phil_rep = (M1-M2) ./ sdpooled;

% Siemens / Philips
sdpooled = sqrt((SD1^2 + SD3^2) / 2);
d_Siem_GE_rep = (M1-M3) ./ sdpooled;

% Siemens / Philips
sdpooled = sqrt((SD2^2 + SD3^2) / 2);
d_Phil_GE_rep = (M2-M3) ./ sdpooled;

%% Demographic table 
DiscoveryTable = MainTable_g1(MainTable_g1.FD < 0.2 & MainTable_g1.FrameTotal >=600 & MainTable_g1.MissingRestSubjects == 0,:);
ReplicationTable = MainTable_g2(MainTable_g2.FD < 0.2 & MainTable_g2.FrameTotal >=600 & MainTable_g2.MissingRestSubjects == 0,:);

Females_d = length(find(DiscoveryTable.demo_gender_id_v2_l_ToolboxTable_g1_discovery == 2));
Males_d = length(find(DiscoveryTable.demo_gender_id_v2_l_ToolboxTable_g1_discovery == 1));
Females_r = length(find(ReplicationTable.demo_gender_id_v2_l_ToolboxTable_g2_replication == 2));
Males_r = length(find(ReplicationTable.demo_gender_id_v2_l_ToolboxTable_g2_replication == 1));

Age_d = mean(DiscoveryTable.tlfb_age_calc_inmonths_l_CBCLTable_g1_discovery);
stdage_d = std(DiscoveryTable.tlfb_age_calc_inmonths_l_CBCLTable_g1_discovery);
min_age_d = min(DiscoveryTable.tlfb_age_calc_inmonths_l_CBCLTable_g1_discovery);
max_age_d = max(DiscoveryTable.tlfb_age_calc_inmonths_l_CBCLTable_g1_discovery);

Age_r = mean(ReplicationTable.tlfb_age_calc_inmonths_l_CBCLTable_g2_replication);
stdage_r = std(ReplicationTable.tlfb_age_calc_inmonths_l_CBCLTable_g2_replication);
min_age_r = min(ReplicationTable.tlfb_age_calc_inmonths_l_CBCLTable_g2_replication);
max_age_r = max(ReplicationTable.tlfb_age_calc_inmonths_l_CBCLTable_g2_replication);


FD_d = mean(DiscoveryTable.FD);
stdFD_d = std(DiscoveryTable.FD);
minFD_d = min(DiscoveryTable.FD);
maxFD_d = max(DiscoveryTable.FD);
FD_r = mean(ReplicationTable.FD);
stdFD_r = std(ReplicationTable.FD);
minFD_r = min(ReplicationTable.FD);
maxFD_r = max(ReplicationTable.FD);

Frametotal_d = mean(DiscoveryTable.FrameTotal);
stdFrametotal_d = std(DiscoveryTable.FrameTotal);
minFrametotal_d = min(DiscoveryTable.FrameTotal);
maxFrametotal_d = max(DiscoveryTable.FrameTotal);
Frametotal_r = mean(ReplicationTable.FrameTotal);
stdFrametotal_r = std(ReplicationTable.FrameTotal);
minFrametotal_r = min(ReplicationTable.FrameTotal);
maxFrametotal_r = max(ReplicationTable.FrameTotal);

NIHtlbxtotal_d = nanmean(DiscoveryTable.nihtbx_totalcomp_uncorrected);
stdNIHtlbxtotal_d = nanstd(DiscoveryTable.nihtbx_totalcomp_uncorrected);
NIHtlbxtotal_r = nanmean(ReplicationTable.nihtbx_totalcomp_uncorrected);
stdNIHtlbxtotal_r = nanstd(ReplicationTable.nihtbx_totalcomp_uncorrected);




