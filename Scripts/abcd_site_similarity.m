function [similarity_matrix, partitions, site_order] = abcd_site_similarity(rmats,site)
% Rum similarity analysis across all abcd subjects
%   Input: rmats = subject correlation matrices (ROI,ROI,subject)
%          site = numeric vector containing site IDs (1-21)

Corrs_by_site_mat = zeros(size(rmats,1)^2/2 - size(rmats,1)/2,size(rmats,3));
idx = find(triu(rmats(:,:,1),1));
for s = 1:size(rmats,3)
    Thisrmat = squeeze(rmats(:,:,s));
    Corrs_by_site_mat(:,s) = Thisrmat(idx);
end

similarity_matrix = FisherTransform(corr(Corrs_by_site_mat));

[a, site_order] = sort(site);

% Find site partitions
SiteDivide = zeros(length(site_order),1);
for r = 1:length(site_order)-1
    if site_order(r) > site_order(r+1) % then new network
        SiteDivide(r,1) = 1;
    end
end
partitions = find(SiteDivide==1);

