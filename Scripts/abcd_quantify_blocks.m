function MeanMatrix = abcd_quantify_blocks(matrix,partitions,roi_order)
% Calculates a mean within blocks of a matrix
%   Input: matrix = an n x n matrix; zeros or ones on the diagonal
%          partitions = vector containing start index of each block; should NOT have a 1 for first block 
%          roi_order = resorted ROIs, to match partitions 

if partitions(1) == 1
        error('Remove the 1 in your partition index!')
end

temp = ones(length(partitions)+1,1);
temp(2:end) = partitions;
partitions = temp;
partitions(end+1) = size(matrix,1);

% Reorder matrix to align with partitions
matrix = matrix(roi_order,roi_order);

% make diagonal nans 
val = diag(matrix);
matrix(matrix == val) = nan;
MeanMatrix = zeros(length(partitions)-1);
for x = 1:length(partitions)-1
    % on diagonal blocks
    Thisstartidx = partitions(x);
    Thisendidx = partitions(x+1)-1;
    % mean 
    MeanMatrix(x,x) = nanmean(nanmean(matrix(Thisstartidx:Thisendidx,Thisstartidx:Thisendidx)));
    
    % off diagonal blocks 
    for y = x+1:length(partitions)-1
        Thisstartidx_x = partitions(x);
        Thisendidx_x = partitions(x+1)-1;
        Thisstartidx_y = partitions(y);
        Thisendidx_y = partitions(y+1)-1;
        MeanMatrix(x,y) = mean(nanmean(matrix(Thisstartidx_x:Thisendidx_x,Thisstartidx_y:Thisendidx_y)));
    end
    diagonal = diag(MeanMatrix);
end
diagonal = diag(MeanMatrix);
MeanMatrix = MeanMatrix + MeanMatrix';
for i = 1:size(MeanMatrix)
   MeanMatrix(i,i) = MeanMatrix(i,i) - diagonal(i);
end

