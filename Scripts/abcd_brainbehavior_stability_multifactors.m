%% All ABCD brain-behavior correlations at the network level 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/ABCD/Scripts/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/scripts/'))
addpath(genpath('/data/nil-bluearc/GMT/Scott/NLAtoolbox/'))
%% Create the data table containing sub ids, corr mats, and all behvioral vars 
abcd_make_big_data_table

%% Extract correlation matrices
% Correlation matrices from the discovery and replication sets  
rmats = cell2mat(MainTable_g1.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g1,1),333);
rmats_disc = permute(rmats,[1,3,2]);

rmats = cell2mat(MainTable_g2.Corrmats);
rmats = reshape(rmats,333,size(MainTable_g2,1),333);
rmats_rep = permute(rmats,[1,3,2]);



%% Get all behavioral variables from both datasets and concatenate

% 1. FD
factors_disc = MainTable_g1.FD;
factors_rep = MainTable_g2.FD;
% 2. Picture vocab
factors_disc(:,2) = MainTable_g1.nihtbx_picvocab_uncorrected;
factors_rep(:,2) = MainTable_g2.nihtbx_picvocab_uncorrected;
% 3. Flanker
factors_disc(:,3) = MainTable_g1.nihtbx_flanker_uncorrected;
factors_rep(:,3) = MainTable_g2.nihtbx_flanker_uncorrected;
% 4. List 
factors_disc(:,4) = MainTable_g1.nihtbx_list_uncorrected;
factors_rep(:,4) = MainTable_g2.nihtbx_list_uncorrected;
% 5. Card sort 
factors_disc(:,5) = MainTable_g1.nihtbx_cardsort_uncorrected;
factors_rep(:,5) = MainTable_g2.nihtbx_cardsort_uncorrected;
% 6. Pattern
factors_disc(:,6) = MainTable_g1.nihtbx_pattern_uncorrected;
factors_rep(:,6) = MainTable_g2.nihtbx_pattern_uncorrected;
% 7. Picture
factors_disc(:,7) = MainTable_g1.nihtbx_picture_uncorrected;
factors_rep(:,7) = MainTable_g2.nihtbx_picture_uncorrected;
% 8. Reading 
factors_disc(:,8) = MainTable_g1.nihtbx_reading_uncorrected;
factors_rep(:,8) = MainTable_g2.nihtbx_reading_uncorrected;
% 9. Fluid Intell
factors_disc(:,9) = MainTable_g1.nihtbx_fluidcomp_uncorrected;
factors_rep(:,9) = MainTable_g2.nihtbx_fluidcomp_uncorrected;
% 10. Crystal Intell 
factors_disc(:,10) = MainTable_g1.nihtbx_cryst_uncorrected;
factors_rep(:,10) = MainTable_g2.nihtbx_cryst_uncorrected;
% 11. Total toolbox
factors_disc(:,11) = MainTable_g1.nihtbx_totalcomp_uncorrected;
factors_rep(:,11) = MainTable_g2.nihtbx_totalcomp_uncorrected;
% 12. Anx/dep
factors_disc(:,12) = MainTable_g1.cbcl_scr_syn_anxdep_r;
factors_rep(:,12) = MainTable_g2.cbcl_scr_syn_anxdep_r;
% 13. Somatic 
factors_disc(:,13) = MainTable_g1.cbcl_scr_syn_somatic_r;
factors_rep(:,13) = MainTable_g2.cbcl_scr_syn_somatic_r;
% 14. Social
factors_disc(:,14) = MainTable_g1.cbcl_scr_syn_social_r;
factors_rep(:,14) = MainTable_g2.cbcl_scr_syn_social_r;
% 15. Thought 
factors_disc(:,15) = MainTable_g1.cbcl_scr_syn_thought_r;
factors_rep(:,15) = MainTable_g2.cbcl_scr_syn_thought_r;
% 16. Attention 
factors_disc(:,16) = MainTable_g1.cbcl_scr_syn_attention_r;
factors_rep(:,16) = MainTable_g2.cbcl_scr_syn_attention_r;
% 17. Rule breaking
factors_disc(:,17) = MainTable_g1.cbcl_scr_syn_rulebreak_r;
factors_rep(:,17) = MainTable_g2.cbcl_scr_syn_rulebreak_r;
% 18. Aggresive 
factors_disc(:,18) = MainTable_g1.cbcl_scr_syn_aggressive_r;
factors_rep(:,18) = MainTable_g2.cbcl_scr_syn_aggressive_r;
% 19. Internalizing 
factors_disc(:,19) = MainTable_g1.cbcl_scr_syn_internal_r;
factors_rep(:,19) = MainTable_g2.cbcl_scr_syn_internal_r;
% 20. Externaling 
factors_disc(:,20) = MainTable_g1.cbcl_scr_syn_external_r;
factors_rep(:,20) = MainTable_g2.cbcl_scr_syn_external_r;
% 21. Total CBCL
factors_disc(:,21) = MainTable_g1.cbcl_scr_syn_totprob_r;
factors_rep(:,21) = MainTable_g2.cbcl_scr_syn_totprob_r;



%% Concatenate brain data and behavioral variables across datasets 

% Filter out high motion subs, missing data subs, & 8 mins of rest data

% discovery set 
rmats_disc = rmats_disc(:,:,MainTable_g1.FD < 0.20 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingBehavioral_g1 == 0 & MainTable_g1.MissingRestSubjects == 0);
factors_disc = factors_disc(MainTable_g1.FD < 0.20 & MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingBehavioral_g1 == 0 & MainTable_g1.MissingRestSubjects == 0,:);
% replication set
rmats_rep = rmats_rep(:,:,MainTable_g2.FD < 0.20 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingBehavioral_g2 == 0 & MainTable_g2.MissingRestSubjects == 0);
factors_rep = factors_rep(MainTable_g2.FD < 0.20 & MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingBehavioral_g2 == 0 & MainTable_g2.MissingRestSubjects == 0,:);


% Concatenated disc and rep corrmats and behavior variables 
allmats = cat(3,rmats_disc,rmats_rep);
allfactors = [factors_disc;factors_rep];
% Reorder matrices according to Gordon/Laumann network partitions
reorder_gordon_laumann_parcels
% Get ROI order 
roi_order = NetworksOrdered(:,1);
%Reorder
allmats_ordered = allmats(roi_order,roi_order,:);

%% Run brain-behavior correlations for each network pair across all factors 
binsize = 25; % bin width of subjects to loop thru
iter = 1000; % number of loops at each bin width
NetworkCorrs = zeros(length(partitionidx)+1,length(partitionidx)+1,floor(size(allfactors,1)/binsize),iter,size(allfactors,2));
NetworkPvals = zeros(length(partitionidx)+1,length(partitionidx)+1,floor(size(allfactors,1)/binsize),iter,size(allfactors,2));
for x = 1:(length(partitionidx)+1)
    for y = x:(length(partitionidx)+1)
        disp(['On network pair ' num2str(x) ' , ' num2str(y)])
        [TheseNetworkCorrs,TheseNetworkPvals] = abcd_blockwise_correlation_iterative_reliability(allmats_ordered,allfactors,partitionidx,binsize,iter,x,y);
        NetworkCorrs(x,y,:,:,:) = TheseNetworkCorrs;
        NetworkPvals(x,y,:,:,:) = TheseNetworkPvals;
    end
    save('/data/nil-bluearc/GMT/Scott/ABCD/NetworkCorrs_21vars_bynetwork_3ksubs.mat','NetworkCorrs','-v7.3');
    save('/data/nil-bluearc/GMT/Scott/ABCD/NetworkPvals_21vars_bynetwork_3ksubs.mat','NetworkPvals','-v7.3');
    clear TheseNetworkCorrs TheseNetworkPvals
end
% Save output
%save('/data/nil-bluearc/GMT/Scott/ABCD/NetworkCorrs_21vars_bynetwork_3ksubs.mat','NetworkCorrs','-v7.3');
%save('/data/nil-bluearc/GMT/Scott/ABCD/NetworkPvals_21vars_bynetwork_3ksubs.mat','NetworkPvals','-v7.3');



