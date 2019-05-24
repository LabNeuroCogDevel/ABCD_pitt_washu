%% Plotting 
% Reorder ROIs from gordon laumann parcellation & create partition index
reorder_gordon_laumann_parcels
% Get ROI order to pass into plotting function 
roi_order = NetworksOrdered(:,1);
%Plot it
plot_adj_matrix(matrix,roi_order,partitionidx);

caxis([-.4 1]);
colormap(jet);