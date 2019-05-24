% This script will generate networks ordered according to Gordon et al 2017
% Output: Partition index = where to draw lines in adj matrix
%         Networks Ordered = ROI indices to order network. Ordered by X
%         coordinate within networks 

addpath(genpath('/data/nil-bluearc/GMT/Scott/ABCD_DCN_CognitionPaper/Scripts'))
NetworksOrdered = [];
Net_labels={'DM';'Vis';'FP';'DA';'VA';'Sal';'CO';'SMH';'SMM';'AUD';'ParMem';'Context';'NONE'};
% Read in parcel IDs, sort, and then read coordinates 
Read_Gordon_Parcel_IDs
[~, NetworksOrdered]  = sort(NetworkIds(:,1));
Read_Gordon_Coordinates


% Find network partitions
NetDivide = zeros(length(NetworksOrdered),1);
for r = 1:length(NetworksOrdered)-1
    if NetworksOrdered(r) > NetworksOrdered(r+1) % then new network
        NetDivide(r,1) = 1;
    end
end
PartitionIdx = find(NetDivide==1);

% Append x and y coords  
for r = 1:length(NetworksOrdered(:,1))
   ThisROI = NetworksOrdered(r,1);
   NetworksOrdered(r,2) = Xcoord(ThisROI);
end

%Order within network
for p = 1:length(PartitionIdx)
    if p == 1
        Thispartition = PartitionIdx(p);
        [Sortedcoords,SortedcoordsIdx] = sort(NetworksOrdered(1:Thispartition,2));
        NewOrder = NetworksOrdered(1:Thispartition,1);
        NewOrder = NewOrder(SortedcoordsIdx);
        NetworksOrdered(1:Thispartition,1) = NewOrder;
        NetworksOrdered(1:Thispartition,2) = Sortedcoords;  
    else
        Thispartition = PartitionIdx(p);
        [Sortedcoords,SortedcoordsIdx] = sort(NetworksOrdered(PartitionIdx(p-1)+1:Thispartition,2));
        NewOrder = NetworksOrdered(PartitionIdx(p-1)+1:Thispartition,1);
        NewOrder = NewOrder(SortedcoordsIdx);
        NetworksOrdered(PartitionIdx(p-1)+1:Thispartition,1) = NewOrder;
        NetworksOrdered(PartitionIdx(p-1)+1:Thispartition,2) = Sortedcoords;  
    end
end