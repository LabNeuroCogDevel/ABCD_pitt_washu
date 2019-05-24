function [ReliabilityCoeffs_edgewise, ReliabilityCoeffs_network, ReliabilityCoeffs_blocks] = abcd_xsub_brainbehav_reliability(FHBrBx,corrmats,behavior,bin,iter,partitions,roi_order)
% Compute  cross subject brain behavior reliability 
% Input: FHBrBx = Replication set brain behavior correlation matrix; should be square  
%        corrmats = correlation matrices (ROI,ROI,subject)
%        behavior = vector of values for behavior (1 per subejct)
%        bin = loop thru subjects in intervals of bin (e.g., 25)
%        iter = number of iterations to calculate split half reliability 
%               at each quantity of subjects
%        partitions = indices where each network begins. DO NOT INCLUDE A 1
%        roi_order = rois ordered by network 

% Get total number of subjects
subtotal = size(corrmats,3);

% Get blockwise mean of first half of data to compare later
MeanMatrix_FH = abcd_quantify_blocks(FHBrBx,partitions,roi_order);

% new partitions
newpartitions= 1;
newpartitions(2:length(partitions)+1,1) = partitions;
newpartitions(length(newpartitions)+1) = size(corrmats,1);

% Get blockwise vectors of first half of data (FHvals) to compare later
for x = 1:length(newpartitions)-1
    % on diagonal blocks
    Thisstartidx = newpartitions(x);
    Thisendidx = newpartitions(x+1)-1;
    tmpidx = find(triu(FHBrBx(Thisstartidx:Thisendidx,Thisstartidx:Thisendidx),1));
    tmpvals = FHBrBx(tmpidx);
    FHvals{x,x,:} = tmpvals;
    % off diagonal blocks 
    for y = x+1:length(newpartitions)-1
        Thisstartidx_x = newpartitions(x);
        Thisendidx_x = newpartitions(x+1)-1;
        Thisstartidx_y = newpartitions(y);
        Thisendidx_y = newpartitions(y+1)-1;
        tmpmat = FHBrBx(Thisstartidx_x:Thisendidx_x,Thisstartidx_y:Thisendidx_y);
        tmpvals = tmpmat(:);
        FHvals{x,y,:} = tmpvals;
    end
end

% Get upper triangle (between network effects) & diagonal (within network effects)
uidx = find(triu(FHBrBx,1));
netuidx = find(triu(MeanMatrix_FH));


% Initiate output
endidx = size(corrmats,3)-round(size(corrmats,3)-50,-2);
ReliabilityCoeffs_edgewise = zeros(size([bin:bin:(size(corrmats,3)-endidx)]',1),iter);
ReliabilityCoeffs_network = zeros(size([bin:bin:(size(corrmats,3)-endidx)]',1),iter);
ReliabilityCoeffs_blocks = zeros(13,13,size([bin:bin:(size(corrmats,3)-endidx)]',1),iter);

% Loop thru sub sizes in increments of 'bin' 'iter' times 

for subnum = [bin:bin:(size(corrmats,3)-endidx)]
    disp(['On sub size ' num2str(subnum)])
    for i = 1:iter
        % permute sib indices 
        p = randperm(size(corrmats,3));
        idx = p(1:subnum);
        % Half Two
        BrBx_SH = zeros(333);
        TheseMats = corrmats(:,:,idx);
        TheseBehav = behavior(idx);
        % Run edge-wise correlations between brain and behavior across subs
        for j = 1:size(TheseMats,1)
            for k = 1+j:size(TheseMats,2)
                BrBx_SH(j,k) = corr(squeeze(TheseMats(j,k,:)),TheseBehav);
            end
        end
        BrBx_SH = BrBx_SH + BrBx_SH';
        
        % Edgewise reliability
        ReliabilityCoeffs_edgewise(subnum/bin,i) = FisherTransform(corr(FHBrBx(uidx),BrBx_SH(uidx)));
        
        % Network (block-wise) reliability
        
        MeanMatrix_SH = abcd_quantify_blocks(BrBx_SH,partitions,roi_order);
        ReliabilityCoeffs_network(subnum/bin,i) = FisherTransform(corr(MeanMatrix_FH(netuidx),MeanMatrix_SH(netuidx)));
       
        % Reliability across subjects within each block 
        for x = 1:length(newpartitions)-1
            % on diagonal blocks
            Thisstartidx = newpartitions(x);
            Thisendidx = newpartitions(x+1)-1;
            tmpidx = find(triu(BrBx_SH(Thisstartidx:Thisendidx,Thisstartidx:Thisendidx),1));
            vals = BrBx_SH(tmpidx);
            ReliabilityCoeffs_blocks(x,x,subnum/bin,i) = corr(vals,FHvals{x,x,:});
            % off diagonal blocks 
            for y = x+1:length(newpartitions)-1
                Thisstartidx_x = newpartitions(x);
                Thisendidx_x = newpartitions(x+1)-1;
                Thisstartidx_y = newpartitions(y);
                Thisendidx_y = newpartitions(y+1)-1;
                tmpmat = BrBx_SH(Thisstartidx_x:Thisendidx_x,Thisstartidx_y:Thisendidx_y);
                vals = tmpmat(:);
                ReliabilityCoeffs_blocks(x,y,subnum/bin,i) = corr(vals,FHvals{x,y,:});
            end
        end
    end
        
    % Get means and display 
    ThisMean = mean(ReliabilityCoeffs_edgewise(subnum/bin,:));
    disp(['Edgewise Reliability for ' num2str(subnum) ' is ' num2str(ThisMean)])
    ThisMean = mean(ReliabilityCoeffs_network(subnum/bin,:));
    disp(['Network Reliability for ' num2str(subnum) ' is ' num2str(ThisMean)])
    clear BrBX
    
end
    
end
        




