function [NetworkCorrs,NetworkPvals] = abcd_blockwise_correlation_iterative_reliability(rmats,factors,partitions,binsize,iter,netidone,netidtwo)
%UNTITLED Summary of this function goes here
%   Input:  rmats = roi x roi x subject matrix. ROIs MUST BE ORDERED BY
%           NETWORK FOLLOWING YOUR PARTITIONS VECTOR!!
%           factors = behavioral item(s) to correlate with RSFC
%                     Should be n x f, where n = subject length and 
%                     f = behavioral variables
%           partitions = indicies marking end of network - should not
%           include a 1 or last item in ROI set (e.g., 333 for Gordon
%           parcellation). 
%           binsize = size of subject bins (intervals of x subjects)
%           iter = number of iterations at each bin size
%           netidone = ID of first network (e.g., 1 = Default)
%           netidtwo = ID of second network (e.g., 4 = Dorsal Attn)
%   Output: NetworkCorrs = subject bin x iteration. Each entry is an
%           RSFC/behavior correlation. 
%           NetworkPvals = subject bin x iteration. Each entry is a pvalue
%           for the associated RSFC/behavior correlation in NetworkCorrs.

subnum = size(rmats,3);
varnum = size(factors,2);

% Create new partitions that include a 1 and last ROI (e.g., 333). 
if partitions(1) == 1
        error('Remove the 1 in your partition index!')
end

temp = ones(length(partitions)+1,1);
temp(2:end) = partitions;
partitions = temp;
partitions(end+1) = size(rmats,1);

% Obtain ROI indices of network 1 and network 2.
First_network_rois = [partitions(netidone):partitions(netidone+1)-1];
Second_network_rois = [partitions(netidtwo):partitions(netidtwo+1)-1];

% Index rmats for these networks
sub_rmats = rmats(First_network_rois,Second_network_rois,:);


% Vectorize correlations (results in correlations x subject)
if size(sub_rmats,1) == size(sub_rmats,2) % within network 
    uidx = find(triu(sub_rmats(:,:,1),1)); % grab upper tri only
    Vecmats = zeros(length(uidx),size(sub_rmats,3));
    for s = 1:size(sub_rmats,3)
        Thismat = sub_rmats(:,:,s);
        Vecmats(:,s) = Thismat(uidx); % create ROIpair x subj matrix
    end
else % between network 
    Vecmats = reshape(sub_rmats,size(sub_rmats,1).*size(sub_rmats,2),size(sub_rmats,3)); 
end

% Loop thru subjects in interval of bin size, iter times at each bin.
NetworkCorrs = zeros(length([binsize:binsize:subnum]),iter,varnum);
NetworkPvals = zeros(length([binsize:binsize:subnum]),iter,varnum);
for i = [binsize:binsize:subnum]
    disp(['On binsize: ' num2str(i)])
    for j = 1:1:iter
        %generate a random vector of length subject
        p = randperm(size(Vecmats,2));
        % grab the first binned size number of subjects
        idx = p(1:i);
        % Index out matrices and factor
        TheseMats = Vecmats(:,idx);
        for f = 1:varnum % Loop thru factors
            TheseBehav = factors(idx,f);
            Thismean = mean(TheseMats)';
            [NetworkCorrs(i/binsize,j,f), NetworkPvals(i/binsize,j,f)] = corr(Thismean,TheseBehav);
        end
    end  
%     for f = 1:varnum
%     disp(['Variable ' num2str(f) ' Min = ' num2str(min(NetworkCorrs(i/binsize,:,f)))])
%     disp(['Variable ' num2str(f) ' Max = ' num2str(max(NetworkCorrs(i/binsize,:,f)))])
%     end
%     disp('******************************')
end
end

