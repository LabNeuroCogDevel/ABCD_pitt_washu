function  plot_adj_matrix(matrix,roi_order,partitionidx)
% plot adjacency matrix and draw network lines
% Input:    matrix = adjacency matrix 
%           roi_order = vector containing ROI order

figure;
imagesc(matrix(roi_order,roi_order));colorbar
hold on
%draw lines
for n = 1:length(partitionidx)
    line([1 length(roi_order)],[partitionidx(n) partitionidx(n)],'Color','k','LineWidth',2.5)
    line([partitionidx(n) partitionidx(n)],[1 length(roi_order)],'Color','k','LineWidth',2.5)
end