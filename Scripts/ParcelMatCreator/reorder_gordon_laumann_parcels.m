% Script to generate network order (and sort by x coord) and partition
% index for line drawing on correlation matrices 

addpath(genpath('/data/nil-bluearc/GMT/Scott/ABCD_Brain_Cog_Paper/Scripts'))
NetworksOrdered = [];
Net_labels={'DM';'Vis';'FPN';'DAN';'VAN';'Sal';'CON';'SMH';'SMM';'AUD';'ParMem';'Context';'NONE'};
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
partitionidx = find(NetDivide==1);

% Append x and y coords  
for r = 1:length(NetworksOrdered(:,1))
   ThisROI = NetworksOrdered(r,1);
   NetworksOrdered(r,2) = Xcoord(ThisROI);
   NetworksOrdered(r,3) = Ycoord(ThisROI);
end

%Order within network
for p = 1:length(partitionidx)
    if p == 1
        Thispartition = partitionidx(p);
        [Sortedcoords,SortedcoordsIdx] = sort(NetworksOrdered(1:Thispartition,2));
        NewOrder = NetworksOrdered(1:Thispartition,1);
        NewOrder = NewOrder(SortedcoordsIdx);
        NetworksOrdered(1:Thispartition,1) = NewOrder;
        NetworksOrdered(1:Thispartition,2) = Sortedcoords;  
    else
        Thispartition = partitionidx(p);
        [Sortedcoords,SortedcoordsIdx] = sort(NetworksOrdered(partitionidx(p-1)+1:Thispartition,2));
        NewOrder = NetworksOrdered(partitionidx(p-1)+1:Thispartition,1);
        NewOrder = NewOrder(SortedcoordsIdx);
        NetworksOrdered(partitionidx(p-1)+1:Thispartition,1) = NewOrder;
        NetworksOrdered(partitionidx(p-1)+1:Thispartition,2) = Sortedcoords;  
    end
end